# Core Genon facility
from __future__ import print_function, division
import numpy as np
import math
from scipy.stats import binom
import json
import sys
import os
import pandas as pd
from collections import Counter, defaultdict
import fisher
import scipy.stats as stats # this is for fisher test. In fact, fisher module is much faster than this,
# but `pip install fisher` has not been successful
import pymongo
import json
import itertools
import pickle
import h5py
import tabix
import re
import functools
import utils
#import itertools
#insample_cutoff = 3
np.seterr(divide='raise')

class partialmethod(functools.partial):
    '''
    python2's functools does not have partialmethod
    '''
    def __get__(self, instance, owner):
        if instance is None:
            return self
        return functools.partial(self.func, instance,
                       *(self.args or ()), **(self.keywords or {}))

def partialclass(cls, *args, **kwds):
    '''
    this is used to subclass with some useful defaults
    '''
    class Gene(cls):
        __init__ = partialmethod(cls.__init__, *args, **kwds)

    return Gene

class GeneBase:
    '''
    gene to analyse, including ground truth HPO and inheritance mode
    '''
    def __init__(self, symbol, mode, hpos, **kwargs):
        if type(hpos) is not list:
            msg = "Gene's hpos has to be a list"
            raise AttributeError(msg)
        self.symbol = symbol
        self.mode = mode
        self.hpos = hpos
        for k,v in kwargs.items():
            setattr(self, k, v)

        # some file locations
        self.gene_file = '../data/private/cutoff/raw_{}_gene_hpo.json'.format(symbol)
        self.vcf_file = '../data/public/phasing/vcf/{}_cleaned.vcf'.format(symbol)
        self.vep_vcf_file = '../data/public/phasing/vcf/vep/{}_vep.vcf'.format(symbol)
        self.phase_file = '../data/public/phasing/{}_combined_phase.pkl'.format(symbol)

    @property
    def vcf(self):
        '''
        get vcf coverage
        '''
        if getattr(self, '_vcf', None) is None:
            vcf = utils.get_cover_from_vcf(self.vcf_file)
            self._vcf = vcf
        return self._vcf

    @property
    def phase(self):
        '''
        get phase data
        '''
        if getattr(self, '_phase', None) is None:
            with open(self.phase_file, 'rb') as inf:
                phase = pickle.load(inf)
            self._phase = phase
        return self._phase

    def get_batch_artefacts(self, gp):
        # using binom, and gnomad_af as p, to produce probability to help identify batch specific artefacts
        # lower_bound is there to remove cohorts where there is just one patient
        # zero_gnomad_c_cutoff allows max internal count when gnomad_af is 0
        dt_d = defaultdict(Counter)
        dt_r = defaultdict(Counter)
        cohorts = Counter()
        for k,v in gp['patients'].items():
            cohorts[ self.patient_mini[k]['contact'] ] += 1
            vc = Counter(v['variants'])
            for i in vc:
                dt_d[ self.patient_mini[k]['contact'] ][i] += 1
                if vc[i] > 1:
                    dt_r[ self.patient_mini[k]['contact'] ][i] += 1
        # remove cohorts with count lower than lower_bound
        for k in cohorts.keys():
            if cohorts[k] < self.lower_bound:
                del cohorts[k]
                del dt_d[k]
                if k in dt_r:
                    del dt_r[k]

        # for heterozygous variants
        result_d = defaultdict(list)
        for k1,v1 in dt_d.items():
            n_variants = len(v1)
            for k2,v2 in v1.items():
                if not gp['variants'][k2]['gnomad_af']:
                    if v2 > self.zero_gnomad_c_cutoff:
                        result_d[k1].append(k2)
                    continue
                prob = 1 - binom.cdf(v2-1, cohorts[k1], gp['variants'][k2]['gnomad_af'])
                if prob < self.binom_cutoff / n_variants:
                    #print(k2,prob)
                    result_d[k1].append(k2)
        for k in result_d:
            result_d[k] = set(result_d[k])

        # for homozygous variants
        result_r = defaultdict(list)
        for k1,v1 in dt_r.items():
            n_variants = len(v1)
            for k2,v2 in v1.items():
                if not gp['variants'][k2]['gnomad_hom_af']:
                    if v2 > self.zero_gnomad_c_cutoff:
                        result_r[k1].append(k2)
                    continue
                prob = 1 - binom.cdf(v2-1, cohorts[k1], gp['variants'][k2]['gnomad_hom_af'])
                if prob < self.binom_cutoff / n_variants:
                    #print(k2,prob)
                    result_r[k1].append(k2)  
        for k in result_r:
            result_r[k] = set(result_r[k])
        return {'d':result_d,'r':result_r}

    # remove batch specific artefact variants
    def remove_batch_artefacts(self, gp, bad_vs, mode='all'):
        result = {
                'patients':{},
                'variants':gp['variants'],
                'gene_id':gp['gene_id'],
                'pat_a':gp['pat_a'],
                }
        bad_p = []
        for k1,v1 in gp['patients'].items():
            cohort = self.patient_mini[k1]['contact']
            this_bad_vs = []
            # collect het artefacts
            if mode != 'r' and cohort in bad_vs['d']:
                this_bad_vs += [i for i in v1['variants'] if i in bad_vs['d'][cohort]]
            # collect hom artefacts
            if mode != 'd' and cohort in bad_vs['r']:
                vc = Counter(v1['variants'])
                for k in vc:
                    if vc[k] > 1 and k in bad_vs['r'][cohort]:
                        this_bad_vs.append(k)
            this_vs = [i for i in v1['variants'] if i not in this_bad_vs]
            if this_vs:
                result['patients'][k1] = {
                        'hpo':v1['hpo'],
                        'variants':this_vs,
                        }
        return result
    def get_genotype_phenotype(self):
        '''
        get the relevant genotype phenotype data
        '''
        data = utils.read_files(self.gene_file)
        # remove nc?
        if self.closeness:
            nc_variants = utils.extract_nc_from_vcf(
                    self.vep_vcf_file,
                    self.closeness
                    )
            utils.remove_noncoding(data,nc_variants)
        utils.cleanse_variants(data)
        # remove batch effect?
        if self.binom_cutoff:
            batch_artefacts = self.get_batch_artefacts(
                    data,
                    )
            data  = self.remove_batch_artefacts(data,batch_artefacts)
        return data

    def get_patient_map(self,mode):
        '''
        get what patients in what bins
        note that this mode may be different from self.mode,
        given a chance to do analysis for the 'wrong' mode
        '''
        gn = np.arange(self.grange[0],self.grange[1],self.steps[1])
        ca = np.arange(self.crange[0],self.crange[1],self.steps[0])
        gp = self.get_genotype_phenotype()
        patient_map = defaultdict(list)
        for i in range(len(ca)):
            p = []
            for j in range(len(gn)):
                p,narrow_vs,not_covered_patients = self.get_patients(
                        gp,
                        mode,
                        (gn[j], gn[j]+self.steps[1]),
                        (ca[i], ca[i]+self.steps[0])
                        )
                patient_map[(i,j)] = (p,not_covered_patients)
        return patient_map

    def get_variants(self,vs,mode,gr,cr):
        '''
        Given mode, gr and cr, return variants
        Narrow variants strictly match gr and cr, 
        and can be used for dominant inheritance.
        For recessive, broad_vs also returned, which
        matches the lower bounds for both gr and cr.
        '''
        mode_dict = {'r':'gnomad_hom_af','d':'gnomad_af'}
        narrow_vs = (k for k,v in vs.items() if gr[0]<=v[mode_dict[mode]]<gr[1] and cr[0]<=self.cadd[k]<cr[1])
        broad_vs = tuple()
        if mode == 'r':
            broad_vs = (k for k,v in vs.items() if v[mode_dict[mode]]<gr[1] and self.cadd[k]>=min(cr[0],self.second_cadd_min))
        return (set(narrow_vs),set(broad_vs))

    def get_patients(self, gp, mode, gr, cr):
        '''
        given ranges, find matching patients
        '''
        p = []
        narrow_vs,broad_vs = self.get_variants(gp['variants'],mode,gr,cr)
        if not narrow_vs:
            return ([],narrow_vs,set())

        # get patients not covered on miss_cutoff* narrow_vs
        s = self.vcf.loc[narrow_vs].sum()
        not_covered_patients = set(s[s<len(narrow_vs)*self.miss_cutoff].index)

        if mode == 'd':
            p = [k for k,v in gp['patients'].items() if set(v['variants']) & narrow_vs]
        elif mode == 'r':
            for k,v in gp['patients'].items():
                good = []
                other = []
                for i in v['variants']:
                    if i in narrow_vs:
                        good.append(i)
                    elif i in broad_vs:
                        other.append(i)
                r_bad = [] # for removing linked variants
                if len(good) and len(good+other) > 1:
                    pos = [int(i.split('-')[1]) for i in good+other]
                    if (len(pos) > len(set(pos))) or (max(pos) - min(pos) > self.gap):
                        if self.phase_cutoff:
                            # any of them is hom?
                            if len(set(good)) < len(good):
                                p.append(k)
                                continue
                            # remove hom in other, as hom in other renders the variants in good unnecessary
                            
                            O = [k1 for k1,v1 in Counter(other).items() if v1 == 1]
                            if len(good+O) < 2: continue
                            for i in itertools.combinations(sorted(good+O,key=lambda x: int(x.split('-')[1])),2):
                                if not set(good) & set(i): continue
                                cis_p = self.phase[k].setdefault(i,0)
                                if cis_p < self.phase_cutoff: 
                                    p.append(k)
                                    break
                                    
                        else:
                            p.append(k)
        return (p,narrow_vs,not_covered_patients)
        


class GenonResult:
    '''
    has genon_sums and genon_ratios
    predicted_mode: recessive if positive, dominant if negative.
    The bigger the margin the more confidence in the call
    '''
    def __init__(self):
        self.genon_sum = defaultdict(lambda: defaultdict(dict))
        self.genon_ratio = defaultdict(lambda: defaultdict(dict))
        self.genes = None
        self.predicted_mode = dict()
    def __str__(self):
        s = ''
        for gene in self.genon_sum:
            s += gene
            # write mode
            if self.genes[gene].mode is not None:
                mode = self.genes[gene].mode
                s += ' - Given mode is {}\n'.format(mode)
            else:
                mode_digit = self.predicted_mode[gene]
                if mode_digit > 0:
                    mode = 'r'
                elif mode_digit < 0:
                    mode = 'd'
                else:
                    mode = 'u'
                s += ' - Predicted mode is {}\n'.format(mode)
            # do not want to carry on with mode:u
            if mode == 'u':
                continue
            s += ' ' * (len(gene) + 1)
            # write genon sums with ratios
            if self.genes[gene].hpos is not None:
                s += '- HPOs are given...'
            else:
                s += '- HPOs are predicted...'
            s += '\n'
            for hpo,Sum in sorted(
                    self.genon_sum[gene][mode].items(), 
                    key = lambda x :x[1],
                    reverse = True):
                s += '\t{}\n'.format(hpo)
                s += '\t\tGenon_sum: {}\n'.format(Sum)
                s += '\t\tGenon_ratio: {}\n'.format(
                        self.genon_ratio[gene][mode][hpo]
                        )
        return s
                

class Genon:
    '''
    this is the class to do the major analysis
    '''
    def __init__(self):
        # some file locations
        self.hdf5 = 'patient_map.hdf5'
        self.cadd_file = '../data/public/cutoff/all_cadd.tsv'
        self.patient_mini_file = (
                '../data/private/hpo/patients_hpo_snapshot'
                '_2017-May_mini.tsv'
                )
        self.patient_info_file = (
                '../data/private/hpo/patients_hpo_snapshot'
                '_2017-May.tsv'
                )
        # some parameters
        # do you want to remove non_coding variants?
        # closeness used to remove noncoding variants {closeness} 
        #  away from any exons.
        # if you do not want to remove non_coding variants,
        #  leave it as None
        self.closeness = 4
        # remove batch effect? if yes, set the binom_cutoff. 
        #  if no, leave it as None
        self.binom_cutoff = 1e-10
        # log(0) is meaningless. replace 0 with {replace_zero}
        self.replace_zero = 1e-6
        # If a cohort's size is lower than {lower_bound},
        # do not process the cohort for batch effect
        self.lower_bound = 2
        # zero_gnomad_c_cutoff allows max internal count when gnomad_af is 0
        self.zero_gnomad_c_cutoff = 2
        # if the matching variants are miss-called by more than
        # {miss_cutoff} of all the patients, discard the variants
        self.miss_cutoff = 0.5
        # phase_cutoff used to phase variants
        self.phase_cutoff = .75
        # for recessive, what is the min for the second variant's cadd
        #  when the first variant's cadd is high(er than {second_cadd_min})?
        self.second_cadd_min = 15
        # only analyse hpo with N >= {N}
        # if no phase data is available, one can set gap to 100,
        # to treat close-by variants as cis
        self.gap = 0
        self.N = 100
        # the following hpos will not be analysed
        # they are inheritance modes that we do not know
        self.hpo_mask = ['HP:0000007','HP:0000006','HP:0003745','HP:0000005']
        # steps = (cadd step, gnomad step)
        self.steps = (5, 0.00025)
        # what is the gnomad range for analysis?
        self.grange = (0, 0.01)
        self.crange = (0, 60)
        # if hpos are not given, Genon will try to find out associated hpos
        self.coefficient = 1.

        # mongodb
        self.conn = pymongo.MongoClient()#(host = 'phenotips')
        self.hpo_db = self.conn['hpo']
        self.p_db = self.conn['patients']

        # genes to analyse
        # subclass Gene with useful defaults
        Gene = partialclass(GeneBase, 
                closeness = self.closeness,
                grange = self.grange,
                crange = self.crange,
                cadd = self.cadd,
                steps = self.steps,
                binom_cutoff = self.binom_cutoff,
                lower_bound = self.lower_bound,
                replace_zero = self.replace_zero,
                second_cadd_min = self.second_cadd_min,
                zero_gnomad_c_cutoff = self.zero_gnomad_c_cutoff,
                gap = self.gap,
                miss_cutoff = self.miss_cutoff,
                patient_mini = self.patient_mini,
                phase_cutoff = self.phase_cutoff,
                )
        self.genes = dict(
                ABCA4 = Gene('ABCA4','r',
                    [
                        'HP:0007754',
                        'HP:0000556',
                        'HP:0000505',
                        'HP:0000512',
                        'HP:0000662',
                        'HP:0000510',
                        ]
                    ),
                SCN1A = Gene('SCN1A','d',
                    [
                        'HP:0001250',
                        # HP:0000707 includes dementia HP:0000726, 
                        # which contributes to noise.
                        'HP:0000707',
                        ]
                    ),
                USH2A = Gene('USH2A','r',
                    [
                        'HP:0000510',
                        'HP:0000556',
                        'HP:0000365',
                        'HP:0000505',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                RPGR = Gene('RPGR','r',
                    [
                        'HP:0000548',
                        'HP:0000556',
                        'HP:0000505',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                GUCY2D = Gene('GUCY2D','d',
                    [
                        'HP:0001103',
                        'HP:0000550',
                        'HP:0000548',
                        'HP:0000639',
                        'HP:0000556',
                        'HP:0000505',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                TERT = Gene('TERT','d',
                    [
                        'HP:0005528',
                        'HP:0002754',
                        'HP:0008404',
                        'HP:0000951',
                        'HP:0000234',
                        ]
                    ),
                PROM1 = Gene('PROM1','d',
                    [
                        'HP:0007754',
                        'HP:0000556',
                        'HP:0000505',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                CRB1 = Gene('CRB1','r',
                    [
                        'HP:0000510',
                        'HP:0000556',
                        'HP:0000505',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                CNGB1 = Gene('CNGB1','r',
                    [
                        'HP:0000510',
                        'HP:0000556',
                        'HP:0000505',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                CERKL = Gene('CERKL','r',
                    [
                        'HP:0000510',
                        'HP:0000556',
                        'HP:0000505',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                BBS1 = Gene('BBS1','r',
                    [
                        'HP:0000510',
                        'HP:0000518',
                        'HP:0000556',
                        'HP:0000505',
                        'HP:0000512',
                        ]
                    ),
                IMPG2 = Gene('IMPG2','r',
                    [
                        'HP:0000510',
                        'HP:0000556',
                        'HP:0000505',
                        'HP:0000512',
                        'HP:0000662',
                        ]
                    ),
                )

    # patient info
    @property
    def patient_info(self):
        if getattr(self, '_patient_info', None) is None:
            patient_info = utils.get_snapshot(self.patient_info_file)
            self._patient_info = patient_info
        return self._patient_info

    @property
    def phs(self):
        # get patient numbers for each hpo
        if getattr(self, '_phs', None) is None:
            self._phs = utils.get_phs(self.patient_info)
        return self._phs

    @property
    def patient_mini(self):
        if getattr(self, '_patient_mini', None) is None:
            patient_mini = utils.get_snapshot(self.patient_mini_file)
            # add contact info
            all_p = self.p_db.patients.find(
                    {'external_id': {'$in':patient_mini.keys()}},
                    {'external_id':1,'contact':1}
                    )
            for i in all_p:
                patient_mini[i['external_id']] = {
                        'hpo': patient_mini[i['external_id']],
                        'contact': i['contact']['user_id']
                        }
            self._patient_mini = patient_mini
        return self._patient_mini

    @property
    def cadd(self):
        '''
        get cadd_phred for all variants
        '''
        if getattr(self, '_cadd', None) is None:
            with open(self.cadd_file,'r') as inf:
                cadd = {}
                for row in inf:
                    if row[0] == '#': continue
                    row = row.rstrip().split('\t')
                    v_id = '-'.join(row[:2]+row[2:4])
                    cadd[v_id] = float(row[-1])
            self._cadd = cadd
        return self._cadd
    
    def phenogenon(self,gene,patient_map,hpo):
        # get raw p_a and p_h
        vcf_patients = set(list(gene.vcf))
        raw_p_a = vcf_patients & set(self.patient_info.keys())
        raw_p_h = set([kk for kk,vv in self.patient_info.items() if hpo in vv and kk in vcf_patients])
        # make a np matrix
        shape1,shape2 = [],[]
        for k in patient_map:
            shape1.append(k[0])
            shape2.append(k[1])
        logp_df = np.zeros( (
            len(set(shape1)),
            len(set(shape2))
            ) )
        for k,v in patient_map.items():
            p_a = len(raw_p_a - v[1])
            p_h = len(raw_p_h - v[1])
            p_g = len(v[0])
            if not p_g:
                logp_df[k[0]][k[1]] = 0
                continue
            p_gh = 0
            for p in v[0]:
                if hpo in self.patient_info[p]:
                    p_gh += 1
            pval = fisher.pvalue(
                    p_a - p_h - p_g + p_gh,
                    p_h - p_gh,
                    p_g - p_gh,
                    p_gh
                    ).right_tail
            logp_df[k[0]][k[1]] = -math.log10(pval or 1e-10)
        return logp_df

    def trim_ns(self,ns,ps):
        '''
        trim ns to make sure there's no negative hpo overlaps positive ones
        '''
        result = []
        for hn in ns:
            if hn in self.hpo_mask:
                continue
            bad = 0
            for hp in ps:
                A = hn in [i['id'][0] for i in utils.get_hpo_ancestors(self.hpo_db, hp)]
                B = hp in [i['id'][0] for i in utils.get_hpo_ancestors(self.hpo_db, hn)]
                if A or B:
                    bad = 1
                    break
            if not bad:
                result.append(hn)
        return result

    def get_positive_negative_hpos(self,genon_sums):
        '''
        get positive and negative hpo sets
        '''
        ps,ns = {'r':[],'d':[]},{'r':[],'d':[]}

        # get positive hpos and negative hpos
        for mode in ('r','d'):
            cutf = np.mean(genon_sums[mode].values()) + \
                    self.coefficient * \
                    np.std(genon_sums[mode].values())
            for k,v in genon_sums[mode].items():
                if v > cutf:
                    ps[mode].append(k)
                elif v < cutf:
                    ns[mode].append(k)
            ps[mode] = utils.hpo_minimum_set(self.hpo_db,ps[mode])
            ns[mode] = utils.hpo_minimum_set(self.hpo_db,ns[mode])
            # ns cant be super/subclass of any ps, and can't be in hpo_mask
            ns[mode] = self.trim_ns(ns[mode],ps[mode])
        return ps,ns

    def analyse(self):
        '''
        Genon analyse
        if hpos for a gene is provided, treat the hpos as truely associated
        else, calculate Phenogenon for all hpos, both modes, and find out
        the associated hpos
        '''
        R = GenonResult()
        R.genes = self.genes
        hpos = [i for i,v in self.phs.items() if v >= self.N and i not in self.hpo_mask]
        for gene in self.genes:
            this_hpos = self.genes[gene].hpos or hpos
            for mode in ('r','d'):
                # get patient map
                patient_map = self.genes[gene].get_patient_map(mode)
                for hpo in this_hpos:
                    genon = self.phenogenon(
                            self.genes[gene],
                            patient_map,
                            hpo
                            )
                    genon_sum = sum(genon[:,0])
                    # silence the RuntimeWarning
                    S = np.sum(genon)
                    if S == 0:
                        genon_ratio = 0
                    else:
                        genon_ratio = genon_sum / S

                    # write to result
                    R.genon_sum[gene][mode][hpo] = genon_sum
                    R.genon_ratio[gene][mode][hpo] = genon_ratio

            # are hpos provided? if not, predict
            if not self.genes[gene].hpos:
                ps, ns = self.get_positive_negative_hpos(R.genon_sum[gene])
                # remove ns
                # note that ns is minised, so need to use ps to remove
                # any unwanted hpos
                for mode in ('d','r'):
                    for hpo in hpos:
                        if hpo not in ps[mode]:
                            R.genon_sum[gene][mode].pop(hpo, None)
                            R.genon_ratio[gene][mode].pop(hpo, None)

            # predict inheritance mode
            gxg = {'r':0,'d':0}
            for mode in ('r','d'):
                for hpo in R.genon_sum[gene][mode]:
                    gxg[mode] += R.genon_sum[gene][mode][hpo] * \
                            R.genon_ratio[gene][mode][hpo]

            R.predicted_mode[gene] = gxg['r'] - gxg['d']

        return R

if __name__ == '__main__':
    G = Genon()
    result = G.analyse()
    print(result)
