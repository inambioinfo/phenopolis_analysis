'''
common functions used across hpo analyses
'''
from __future__ import print_function, division
import sys
sys.path.append('../commons')
import phenopolis_utils
import math
from collections import defaultdict
import json
import os

dbs = phenopolis_utils.get_mongo_collections()
'''
information content
'''
def IC(x,t):
    return -math.log(x/t)

'''
get json data
'''
def get_json(input):
    with open(input,'r') as inf:
        return json.load(inf)

'''
get json data
'''
def write_json(data,outfile):
    with open(outfile,'w') as outf:
        json.dump(data, outf)

'''
expand hpo
'''
def expand_hpo(data):
    for d in data:
        all_hpos = [] # union of all ancestors
        for hpo in d['hpo']:
            anc = phenopolis_utils.get_hpo_ancestors(dbs['hpo_db'], hpo)
            for a in anc:
                all_hpos.extend(a['id'])
        d['hpo'] = list(set(all_hpos))
'''
get hpo_freq
'''
def get_hpo_freq(data):
    hpo_freq = defaultdict(int)
    for d in data:
        for h in d['hpo']:
            hpo_freq[h] += 1
    return hpo_freq

'''
ancient wrapper
'''
def get_ancestors(id):
    return [i['id'][0] for i in phenopolis_utils.get_hpo_ancestors(dbs['hpo_db'],id)]

'''
get nearest common ancestor
'''
def get_nearest_common_ancestor(h1,h2,freq):
    hpos = list(set(get_ancestors(h1)) & set(get_ancestors(h2)))
    fs = [freq[h] for h in hpos]
    return hpos[fs.index(min(fs))]
    
'''
parent
'''
def get_parents(id):
    result = []
    h=dbs['hpo_db'].hpo.find_one({'id':id})
    return h.get('is_a',[])

'''
lin's similarity
'''
def lin_similarity(h1,h2,freq):
    anc = get_nearest_common_ancestor(h1,h2,freq)
    t = freq['HP:0000001']
    return 2*IC(freq[anc],t)/( IC(freq[h1],t) + IC(freq[h2],t) )
    
'''
weighted Lin's similarity
'''
def weighted_lin(h1,h2,real_matrix,sim_mean_matrix,freq):
    k = '-'.join(sorted([h1,h2]))
    real_weight = real_matrix.get(k,0)
    sim_weight = sim_mean_matrix.get(k,0)
    return real_weight * lin_similarity(h1,h2,freq)
'''
phenotype similarity between two sets of hpos
method now only supports BM (best match)
'''
def patient_hpo_similarity(hpos1,hpos2,method,kernel):
    # read matrix and freq
    matrix_file = os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['matrix_file'])
    matrix = get_json(matrix_file)
    # get keys and their values in the matrix
    weights = {}
    for h1 in hpos1:
        for h2 in hpos2:
            k = '-'.join(sorted([h1,h2]))
            weights[k] = matri
    # take two directions average
    # hpos1 -> hpos2
    for h1 in hpos1:
        pass
