#!/bin/env python
'''
patient matrix for clustering
'''
from __future__ import print_function, division
import sys
sys.path.append('../commons')
import phenopolis_utils
import hpo_helper
from collections import defaultdict
import os
import itertools

# methods
'''
calulate two sets HPOs' similarity
'''
def similarity(dbs,hpos1,hpos2,freq,matrix,kernel,mode='sum',mx=True):
    # minimise
    hpos1 = phenopolis_utils.hpo_minimum_set(dbs['hpo_db'],hpos1)
    hpos2 = phenopolis_utils.hpo_minimum_set(dbs['hpo_db'],hpos2)
    result1 = result2 = 0
    # h1 ---> h2
    errmsg = '{h} not in freq'
    cache = {}
    for h1 in hpos1:
        if h1 not in freq:
            raise ValueErorr(errmsg.format(h=h1))
        this = []
        for h2 in hpos2:
            if h2 not in freq:
                raise ValueError(errmsg.format(h=h2))
            params = {
                    'h1':h1,
                    'h2':h2,
                    'freq':freq,
                    'matrix':matrix,
            }
            sim = kernel(params)
            cache[(h1,h2)] = sim
            this.append(sim)
        if mx:
            result1 += max(this)
        else:
            result1 += sum(this)

    # h2 ---> h1
    for h2 in hpos2:
        this = []
        for h1 in hpos1:
            this.append(cache[(h1,h2)])
        if mx:
            result2 += max(this)
        else:
            result2 += sum(this)
    
    # average?
    if mode == 'average_cooc':
        result1 = result1 / max(len(hpos1),len(hpos2))
        result2 = result2 / max(len(hpos1),len(hpos2))
    elif mode == 'average_lin':
        result1 = result1 / len(hpos1)
        result2 = result2 / len(hpos2)
    return (result1 + result2) / 2

'''
wrapper of lin_similarity
'''
def lin_similarity(params):
    return hpo_helper.lin_similarity(params['h1'],params['h2'],params['freq'])

'''
using co-occurrence matrix
'''
def cooc_matrix(params):
    k = '-'.join(sorted([params['h1'],params['h2']]))
    return params['matrix'].get(k,0)

'''
marry lin and cooc together
'''
def marry(params):
    return lin_similarity(params) * cooc_matrix(params)
'''
produce a matrix for patients
mode = ['sum','average']
'''
def patient_hpo_matrix(dbs,snapshot,freq,matrix,kernel,mode='sum',mx=True):
    result = {}
    for k in itertools.combinations(snapshot.keys(),2):
        print(k)
        result[tuple(sorted(k))] = similarity(dbs,snapshot[k[0]]['hpos'],snapshot[k[1]]['hpos'],freq,matrix,kernel,mode,mx)
    return result

if __name__ == '__main__':
    # get db and snapshot and freq and matrix
    dbs = phenopolis_utils.get_mongo_collections()
    snapshot = hpo_helper.read_snapshot()
    freq_file = os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['mini_real_freq_file'])
    freq = hpo_helper.get_json(freq_file)

    real_matrix_file = os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['matrix_file'])
    real_matrix = hpo_helper.get_json(real_matrix_file)

    test = {
            'brain_p1': snapshot['Vulliamy_April2014_Sample_3482'],
            'brain_p2': snapshot['Vulliamy_April2014_Sample_3615'],
            'myopa_p1': snapshot['WebsterURMD_Sample_05G06386'],
            'myopa_p2': snapshot['WebsterURMD_Sample_05G03129'],
            'ush_p1': snapshot['IRDC_batch5_OXF_3023'],
            'ush_p2': snapshot['IRDC_batch3_MAN_1020_10003846'],
    }
    h1 = ['HP:0000365','HP:0010978']
    h2 = ['HP:0000365','HP:0010978']
    # cooc
    #matrix = patient_hpo_matrix(dbs,test,freq,real_matrix,cooc_matrix,mode='average')
    sim = similarity(dbs,h1,h2,freq,real_matrix,cooc_matrix,mode='average',mx=True)
    print(sim)
