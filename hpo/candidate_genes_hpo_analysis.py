#!/bin/env python
'''
use candidate genes and solved patients for hpo analysis
'''
from __future__ import print_function, division
import sys
sys.path.append('../commons')
import phenopolis_utils
import hpo_helper
import patient_hpo_matrix
import os
import json
import random

'''
subset genes given cutoff
'''
def subset(genes,cutoff):
    result = {}
    for k,v in genes.items():
        A = len(v['data']) >= cutoff
        B = 1 in [d['solve'] for d in v['data']]
        if A or B:
            result[k] = v
    return result

'''
slim
'''
def slim(genes):
    for k,v in genes.items():
        for d in v['data']:
            if 'genes' in d: 
                del d['genes']
                d['hpo'] = phenopolis_utils.hpo_minimum_set(dbs['hpo_db'],[h['id'] for h in d['hpo'] if h['observed'] == 'yes'])
                d['contact'] = d['contact']['user_id']

'''
fmt
'''
def fmt(genes):
    result = {}
    for k,v in genes.items():
        for d in v['data']:
            result[d['external_id']] = {
                    'gene':k,
                    'hpos':d['hpo'],
                    'symbol':v['symbol'],
                    'contact':d['contact'],
            }
    return result

'''
clean patients
'''
def clean(patients, freq):
    for p in list(patients.keys()):
        if p not in freq:
            del patients[p]

'''
expand
'''
def expand(ps, snapshot, size=None):
    size = size or len(ps)
    # remove patients with HP:0000001 and HP:0000118
    bad_hpos = {'HP:0000001','HP:0000118'}
    bad_keys = [k for k,v in snapshot.items() if bad_hpos & set(v['hpos'])]
    ks = set(snapshot.keys()) - set(ps.keys()) - set(bad_keys)
    for p in random.sample(ks, size):
        # get contact
        contact = dbs['patient_db'].patients.find_one({'external_id':p})['contact']['user_id']
        ps[p] = {
                'gene':'unknown',
                'hpos':snapshot[p]['hpos'],
                'symbol':'unknown',
                'contact':contact,
        }
if __name__ == '__main__':
    dbs = phenopolis_utils.get_mongo_collections()
    genes = phenopolis_utils.get_candidate_genes(dbs)
    # using cutoff's parameter to subset genes
    cutoff = phenopolis_utils.OFFLINE_CONFIG['cutoff']['unsolved_common_min']
    # get freq
    freq_file = os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['mini_real_freq_file'])
    freq = hpo_helper.get_json(freq_file)
    genes = subset(genes,cutoff)
    # slim genes
    slim(genes)
    # some patients are annotated badly. what happens if ...
    # say ABCA4 patients all have 'autosomal recessive..' and 'macular dystrophy' ?
    abca4fake = False
    if abca4fake:
        for g in genes['ENSG00000198691']['data']:
            g['hpo'] = phenopolis_utils.hpo_minimum_set(dbs['hpo_db'],g['hpo'] + ['HP:0000007'])

    # convert genes into patients
    patients = fmt(genes)
    # remove patients if they are not unrelated
    snapshot = hpo_helper.read_snapshot()
    clean(patients,snapshot)
    # randomly choose equal number of unsolved patients who have phenotypes (not just 'all'), add them to patients
    random.seed('abc123')
    expand(patients,snapshot,30)
    f = '../data/private/hpo/patient_info.json'
    hpo_helper.write_json(patients,f)
    # get patient matrix
    real_matrix_file = os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['matrix_file'])
    real_matrix = hpo_helper.get_json(real_matrix_file)
    
    matrix = patient_hpo_matrix.patient_hpo_matrix(dbs,patients,freq,real_matrix,patient_hpo_matrix.cooc_matrix,mode='average_cooc',mx=True)
    # convert keys
    for k in list(matrix.keys()):
        newk = ','.join(k)
        matrix[newk] = matrix.pop(k)
    
    f = '../data/private/hpo/patient_hpo_matrix_cooc_MA.json'
    hpo_helper.write_json(matrix,f)
