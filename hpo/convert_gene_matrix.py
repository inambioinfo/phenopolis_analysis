#!/bin/env python
'''
group patient hpo matrices
'''
from __future__ import print_function, division
import sys
sys.path.append('../commons')
import phenopolis_utils
import hpo_helper
import patient_hpo_matrix
import os
import json
import igraph
import collections

if __name__ == '__main__':
    infile = os.path.join('../data/public/hpo/gene_matrix.json')
    matrix = hpo_helper.get_json(infile)
    genes = set()
    for k in list(matrix.keys()):
        ks = k.split(',')
        genes = genes | set(ks)
        matrix[tuple(ks)] = matrix.pop(k)
    genes = list(genes)
    # build g graph
    g = igraph.Graph([(genes.index(i[0]),genes.index(i[1])) for i in sorted(list(matrix.keys()))])
    g.vs["id"] = genes
    g.es["weight"] = [matrix[i] for i in sorted(list(matrix.keys()))]
    # get community, using fastgreedy or walktrap
    community = g.community_fastgreedy(weights='weight') #can use community_walktrap
    #community = g.community_walktrap(weights='weight')
    group_counts = community.optimal_count
    groups = community.as_clustering(group_counts)
    merges = community.merges
    # make it ready for dendrogram
    nodes = []
    pool = range(len(merges)*4)
    def recur(nodes,j,this_node):
        if j < len(g.vs['id']):
            this_node = '.'.join([this_node,g.vs['id'][j]])
            nodes.append(this_node)
            return nodes
        this_node = '.'.join([this_node,'node'+str(pool.pop(0))])
        nodes.append(this_node)
        for k in merges[j - len(g.vs['id'])]:
            recur(nodes,k,this_node)

    new_node='master_node'
    for v in merges[-1]:
        recur(nodes,v,new_node)

    # write csv
    outfile = '../data/public/hpo/gene_matrix.csv'
    header = ['id','gene']
    with open(outfile,'w') as outf:
        outf.write(','.join(header) + '\n')
        outf.write(','.join(['master_node','']) + '\n')
        for n in nodes:
            row = [n,'']
            end_node = n.split('.')[-1]
            if not end_node.startswith('node'):
                row[1] = ''
            outf.write(','.join(row) + '\n')
    print('done')
