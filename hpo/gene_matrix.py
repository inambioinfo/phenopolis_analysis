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
    result1 = result2 = 0
    # g1 --> g2
    for p1 in g1:
        this = []
        for p2 in g2:
            k = ','.join(sorted([p1,p2]))
            this.append(matrix[k])
        result1 += np.mean(this)

    # g2 --> g1
    for p2 in g2:
        this = []
        for p1 in g1:
            k = ','.join(sorted([p1,p2]))
            this.append(matrix[k])
        result2 += np.mean(this)

    result1 = result1/len(g1)
    result2 = result2/len(g2)
    return (result1 + result2) / 2

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
    matrix = gene_matrix(genes,patient_matrix)
    outfile = '../data/public/hpo/gene_matrix.json'
    hpo_helper.write_json(matrix,outfile)
