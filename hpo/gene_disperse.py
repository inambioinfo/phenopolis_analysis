#!/bin/env python
'''
use some patients with candidate genes to test the performance of MBA-WAM
'''
from __future__ import print_function, division
import sys
sys.path.append('../commons')
import phenopolis_utils
import hpo_helper
import os
import json
import numpy
import pandas as pd
from scipy import stats

'''
get disperse p
'''
def disperse(matrix,patients,g,mode):
    g_patients = [p for p,v in patients.items() if v['symbol'] == g]
    retinal_patients = set([p for p,v in patients.items() if v['contact'] in ['UKIRDC','HARDCASTLE','BLACK','WEBSTER']]) - set(g_patients)
    inn = []
    out = []
    for k,v in matrix.items():
        k = set(k.split(','))
        # not interested in out-of-retinal patients relations if mode == zoom in
        if mode == 'zoom_in' and k - retinal_patients - set(g_patients): continue
        if len(k & set(g_patients)) == 2:
            inn.append({'gene':g,'value':v})
        elif len(k & set(g_patients)) == 1:
            p = list(k - set(g_patients))[0]
            out.append({'gene':patients.get(p,{'symbol':'unknown'})['symbol'],'value':v})
        

    disperse = stats.ttest_ind([i['value'] for i in inn],[i['value'] for i in out],equal_var=False)
    if pd.isnull(disperse[1]):
        return None
    return disperse[1]/2

if __name__ == '__main__':
    #dbs = phenopolis_utils.get_mongo_collections()
    #genes = phenopolis_utils.get_candidate_genes(dbs)

    # get matrix file
    inpath = '../data/private/hpo'
    fs = {'cooc':'patient_hpo_matrix_cooc_MA.json','lin':'patient_hpo_matrix_lin_MA.json'}
    modes = ['zoom_in','zoom_out']
    patients = hpo_helper.get_json('../data/private/hpo/patient_info.json')
    genes = set([v['symbol'] for p,v in patients.items()])
    result = {}
    for fi,fv in fs.items():
        matrix = hpo_helper.get_json(os.path.join(inpath,fv))
        for m in modes:
            result['-'.join([fi,m])] = {}
            for g in genes:
                result['-'.join([fi,m])][g] = disperse(matrix,patients,g,m)
    outfile = '../data/public/hpo/gene_disperse.json'
    hpo_helper.write_json(result,outfile)
