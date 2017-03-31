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

if __name__ == '__main__':
    patient_json = hpo_helper.get_json('../data/private/hpo/patient_info.json')
    # anonymise patients, build a dict
    patients = []
    p_dict = {}
    n = 0
    for k in patient_json:
        p_dict[k] = 'p' + str(n)
        p_dict['p'+str(n)] = k
        patients.append('p'+str(n))
        n += 1

    for f in os.listdir('../data/private/hpo'):
        if not f.startswith('patient_hpo_matrix'): continue

        infile = os.path.join('../data/private/hpo',f)
        matrix = hpo_helper.get_json(infile)
        for k in list(matrix.keys()):
            matrix[tuple(k.split(','))] = matrix.pop(k)
        # build g graph
        g = igraph.Graph([(patients.index(p_dict[i[0]]),patients.index(p_dict[i[1]])) for i in sorted(list(matrix.keys()))])
        g.vs["id"] = patients
        g.vs["gene"] = [patient_json[p_dict[p]]['symbol'] for p in patients]
        g.vs["contact"] = [patient_json[p_dict[p]]['contact'] for p in patients]
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
        outpath = '../data/public/hpo'
        f = f.replace('.json','.csv')
        outfile = os.path.join(outpath,f)
        header = ['id','gene','contact']
        with open(outfile,'w') as outf:
            outf.write(','.join(header) + '\n')
            outf.write(','.join(['master_node','','']) + '\n')
            for n in nodes:
                row = [n,'','']
                end_node = n.split('.')[-1]
                if  end_node[0] == 'p':
                    row[1] = patient_json[p_dict[end_node]]['symbol']
                    row[2] = patient_json[p_dict[end_node]]['contact']
                outf.write(','.join(row) + '\n')
