'''
adapted from phenogenon for a combo of rod-cone and hearing impairment
'''
from __future__ import print_function, division
import tabix
import gnomad_utils
import subprocess
from optparse import OptionParser
import sys
sys.path.append('../commons')
import phenopolis_utils
import os
import json
from collections import defaultdict,Counter
import pandas as pd
import numpy as np
import helper
from sklearn.cluster import KMeans
import itertools
import copy

MONGO = phenopolis_utils.get_mongo_collections()

def get_hpo_from_json(f):
    '''
    if remote server is somehow unavailable, use a local json file instead
    '''
    with open(f,'r') as inf:
        data = '[' + inf.read().rstrip().replace('\n',',') + ']'
        data = json.loads(data)
    # convert it to a dict
    return {i['id'][0]:i for i in data}

def hpo_name(hpo_db, ids):
    '''
    this is only for notebook
    '''
    records = hpo_db.hpo.find({'id':{'$in':ids}},{'id':1,'name':1,'_id':0})
    result = {}
    for r in records:
        result[r['id'][0]] = r['name'][0]
    return result

'''
given chromosome and db, return gene_ranges
'''
def get_chrom_genes(chrom,fields, db):
    # give chrom numbers, get all genes on them

    chrom = str(chrom)
    if chrom not in phenopolis_utils.VALID_CHROMOSOMES:
        raise ValueError('Error: %s is not a valid chromosome!' % chrom)
    gene_ranges = db.genes.find({'chrom':chrom},fields,no_cursor_timeout=True)

    return gene_ranges

'''
when mongodb is not available!
'''
def get_chrom_genes_with_jq(chrom,json_file):
    cmd = """/share/apps/genomics/jq -c '[.gene_id, .gene_name, .chrom, .start, .stop, .xstart, .xstop] | select(.[2]=="%s")|{gene_id:.[0],gene_name:.[1],chrom:.[2],start:.[3],stop:.[4],xstart:.[5],xstop:.[6]}' """ % chrom
    result = subprocess.check_output(cmd+json_file,shell=True)
    return helper.split_iter(result)

def main(**kwargs):
    '''
    parameters:
     genes: optional
     N (selecting HPO with at least N Ph. affecting both \
          #positive (selecting parental HPO in the positive set \
          #and negative set)
     vcf file location
     gnomad files location
     patient_mini, patient_info, both are json files
     cadd path
     unrelated file used to subset vcf file
     v cutoff and p cutoff are to remove variants and patients with \
          #low coverage over the gene

    returns hpo goodness of fit score, p_g (gnomad_freq 
    '''
    # check args
    compulsory_keys = {
        'N',
        'gnomad_path',
        'patient_mini_file',
        'patient_info_file',
        'unrelated_file',
        'gnomad_cutoff',
        'gnomad_step',
        'gnomad_path',
        'hpo_mask',
        'cadd_step',
        'cadd_min',
        'output',
        }
    helper.check_args(compulsory_keys, kwargs, 'main')
    # defaults
    kwargs.setdefault('gene_inheritance_mode',{})
    # output already exist?
    if os.path.isfile(kwargs['output']):
        print('already done')
        return None
    # get patient_mini and patient_info
    patient_info = helper.get_snapshot(kwargs['patient_info_file'])
    patient_mini = helper.get_snapshot(kwargs['patient_mini_file'])
    # get p_h for all hpos
    phs = helper.get_phs(patient_info)
    # add cohort info into patient_mini
    all_p = MONGO['patient_db'].patients.find({'external_id':{'$in':patient_mini.keys()}},{'external_id':1,'contact':1})
    for i in all_p:
        # !!!! this belongs to UCLex's problem!!! remove if publish
        # JingYu and BLACK to UKIRDC, KELSELL to DavidKelsell
        contactdict = dict(
                JingYu = 'UKIRDC',
                Black = 'UKIRDC',
                KELSELL = 'DavidKelsell',
                TonySegal = 'SEGAL',
                SanjaySisodiya = 'SISODIYA',
                )
        contact = i['contact']['user_id']
        contact = contactdict.get(contact,contact)
        patient_mini[i['external_id']] = {'hpo': patient_mini[i['external_id']],
                                          'contact': contact}
    # get hpodb from json
    hpo_db = get_hpo_from_json('../tests/data/new-hpo-hpo.json')
    # get genes, if not provided. get all gene_ids from mongodb, \
            #if provided, convert to gene_id
    fields = {
            'gene_id':1,
            'gene_name':1,
            '_id':0,
            'chrom':1,
            'start':1,
            'stop':1,
            'xstart':1,
            'xstop':1,
            }
    this = {}

    gene_ranges = get_chrom_genes(kwargs['chrom'], fields, MONGO['phenopolis_db'])
    # get gnomad and cadd steps
    gnomad_steps = np.arange(
            0,
            kwargs['gnomad_cutoff']+kwargs['gnomad_step'],
            kwargs['gnomad_step']
            )
    cadd_steps = np.arange(kwargs['cadd_min'], 60, kwargs['cadd_step'])

    # get patient_maps
    with open(
            os.path.join(
                kwargs['patient_maps_path'],
                '{}.json'.format(kwargs['chrom'])
            ),
            'r',
    ) as inf:
        patient_maps = json.load(inf)

    # for each gene, get all valid variants/patients according to p/v_cutoff, 
    # annotate using gnomad
    result = {}
    number_processed = 0
    coding_variants = None
    outf = open(kwargs['output'], 'w')
    outf.write('{')
    for gene_range in gene_ranges:
        # get patient_map
        patient_map = patient_maps.pop(gene_range['gene_id'], None)
        if patient_map is None:
            continue
        # print progress
        number_processed += 1
        if not number_processed % 100:
            print('===processed {} genes==='.format(number_processed))
        print('processing {}'.format(gene_range['gene_name']))

        modes = kwargs['gene_inheritance_mode'].get(
                gene_range['gene_name'],
                kwargs['gene_inheritance_mode'].get(
                    gene_range['gene_id'],
                    'rd'
                    )
                )

        # translate patient_map's key
        pm = {}
        for m in modes:
            pm[m] = {}
            for k,v in patient_map['patient_map'][m].items():
                key = tuple([int(i) for i in k.split(',')])
                pm[m][key] = v
        phenogenon_cache = {'r':{},'d':{}}
        # get phenogenon sums on the first gnomad bin.
        # get all testable hpos
        hpos_to_test = None
        if kwargs['hpos_to_test']:
            hpos_to_test = [i for i in kwargs['hpos_to_test']
                    if i not in kwargs['hpo_mask']]
        else:
            hpos_to_test = [i for i,v in phs.items() 
                    if v >= kwargs['N'] 
                    and i not in kwargs['hpo_mask']]
        for hpos in hpos_to_test:
            # inheritance mode: r and d
            # Note that for each HPO, it only keeps the inheritance mode
            #  with the higher hgf score

            for mode in modes:
                args = dict(
                        hpos = hpos,
                        mode = mode,
                        patient_info = patient_info,
                        patient_map = pm[mode],
                        )

                genon =  helper.phenogenon(**args)

                phenogenon_cache[mode][hpos] = genon.tolist()

        output = json.dumps({
            gene_range['gene_id']:{
                'symbol': gene_range['gene_name'],
                'phenogenon': phenogenon_cache,
                'NP': patient_map['NP'],
            }
        })
        # strip off the braces
        output = output[1:-1]
        if number_processed != 1:
            # meaning not the first record. add a comma
            output = ',' + output
        outf.write(output)

    # close cursor
    gene_ranges.close()

    outf.write('}')

    outf.close()


if __name__ == '__main__':
    # in the end some of the args have to go to the config
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)

    parser.add_option("--chrom",
                      dest="chrom",
                      help="which chrom to process?")
    parser.add_option("--output",
                      dest="output",
                      help="output file name?")
    (options, args) = parser.parse_args()
    args = dict(
        chrom = options.chrom,
        output = options.output,
        gnomad_path = '/cluster/project8/vyp/gnomad_data',
        uclex_genes_json = '../tests/data/uclex-genes.json',
        patient_mini_file = '../data/private/hpo/patients_hpo_'+
            'snapshot_2017-May_mini.tsv',
        patient_info_file = '../data/private/hpo/patients_hpo_'+
            'snapshot_2017-May.tsv',
        patient_maps_path = '../data/private/cutoff/patient_maps_August2017',
        # only look at hpos with patients greater than this number
        N = 60,
        # use this to help find top HPO terms
        # this cutoff is to get poorly covered individuals 
        #  for a given set of variants, to get patient_map
        # e.g. we have a~e five variants. if vmc is set at 0.5,
        #  and an individual not covered on a,b,c, it is removed 
        #  from the analysis

        # this is not implemented yet but essentially adds more weight
        #  to variants with higher cadd scores for recessive mode.
        unrelated_file = '/SAN/vyplab/UCLex/KING/UCL-exome_unrelated.txt',
        # known gene inheritance mode. if provided, no need to infer from data
        #  for sometimes it does make mistakes such that for CERKL
        #gene_inheritance_mode = dict(
        #    ABCA4 = 'r',
        #    CERKL = 'r',
        #    SCN1A = 'd',
        #    GUCY2D = 'd',
        #    USH2A = 'r',
        #    PROM1 = 'd',
        #    TERT = 'd',
        #    CNGB1 = 'r',
        #    CRB1 = 'r',
        #    IMPG2 = 'r',
        #    RPGR = 'r',
        #    ),
        gene_inheritance_mode = {},
        cadd_step = 5,
        cadd_min = 0,
        gnomad_step = 0.00025,
        gnomad_cutoff = 0.01,
        # HPOs not wanted to be included in the analysis
        #  usually the ones people don't record, such as 
        #  inheritance mode
        hpo_mask = (
            'HP:0000007',
            'HP:0000006',
            'HP:0003745',
            'HP:0000005'
            ),
        hpos_to_test = ('HP:0000510,HP:0000365',),
    )
    main(**args)
    print('==done==')
