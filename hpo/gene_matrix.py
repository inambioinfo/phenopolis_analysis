#!/bin/env python
'''
building gene mawtrix from candidate gene patient anlaysis
'''
from __future__ import print_function, division
import sys
sys.path.append('../commons')
import phenopolis_utils
import json
import hpo_helper
from collections import defaultdict
import itertools
import numpy as np

'''
calculate gene-gene similarity
'''
def gene_similarity(g1,g2,matrix):
    result = []
    for p1 in g1:
        for p2 in g2:
            k = ','.join(sorted([p1,p2]))
            result.append(matrix[k])

    return np.mean(result)

'''
group patients into genes
'''
def group(patients):
    result = defaultdict(list)
    for k,v in patients.items():
        result[v['symbol']].append(k)
    return result

'''
build gene_matrix
'''
def gene_matrix(genes, matrix):
    result = defaultdict(float)
    for i in itertools.combinations(genes.keys(),2):
        k = ','.join(sorted(i))
        result[k] += gene_similarity(genes[i[0]],genes[i[1]],matrix)
    return result

if __name__ == '__main__':
    patient_matrix_f = '../data/private/hpo/patient_hpo_matrix_cooc_MA.json'
    patient_matrix = hpo_helper.get_json(patient_matrix_f)
    patients = hpo_helper.get_json('../data/private/hpo/patient_info.json')
    # group them into genes
    genes = group(patients)
    print(gene_similarity(['IRDC_batch2_OXF_3006','IRDC_batch1_LDS_4001_001_543'],['IRDC_batch2_MAN_1007_11002952','Levine_Aug2014_IBDAJ_436','Vulliamy_April2014_Sample_2199'],patient_matrix))
    sys.exit()
    matrix = gene_matrix(genes,patient_matrix)
    outfile = '../data/public/hpo/gene_matrix.json'
    hpo_helper.write_json(matrix,outfile)
