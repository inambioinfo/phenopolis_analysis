#!/bin/env python
'''
generate co-occurrence matrix for HPO [unrelated]
'''
from __future__ import print_function, division
import sys
sys.path.append('../commons')
import phenopolis_utils
import hpo_helper
import math
import json
import os
from collections import defaultdict
import itertools
dbs = phenopolis_utils.get_mongo_collections()

'''
information content
'''
def IC_hpo(id,total,mode,hpo_freq):
    if not id in hpo_freq[mode]:
        return None
    return hpo_helper.IC(len(hpo_freq[mode][id]),total)


'''
snapshot
'''
def get_snapshot():
    with open(os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['snapshot_file']),'r') as inf:
        snapshot = {}
        for row in inf:
            if row[0] == '#': continue
            row = row.rstrip().split('\t')
            snapshot[row[0]] = {'unrelated':int(row[1]),'hpo':row[2].split(',')}
    return snapshot

if __name__ == '__main__':
    # read hpo_freq and snapshot
    hpo_freq = hpo_helper.get_json(os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['hpo_freq_file']))
    snapshot = get_snapshot()

    # get some constant
    total = len(hpo_freq['unrelated']['HP:0000001'])
    IC_max = hpo_helper.IC(1,total)

    # getting serious now!
    result = defaultdict(float)
    #ancestor_cache={}
    method_cache = {}
    mode = 'unrelated'
    n=0
    for k,v in snapshot.items():
        if not v['unrelated']: continue
        n += 1
        print(n)
        print(k)
        hpos = v['hpo']

        # apparently each hpo co-occur with itself
        for h in hpos:
            key = '-'.join([h,h])

            method_cache[h] = method_cache.get(h,IC_hpo(h,total,mode,hpo_freq)/(IC_max))
            result[key] = result.get(key, method_cache[h])
        # check each combination of the hpos 
        for h in itertools.combinations(hpos,2):
            #ancestor_cache[h[0]] = ancestor_cache.get(h[0], hpo_helper.get_ancients(h[0]))
            #ancestor_cache[h[1]] = ancestor_cache.get(h[1], [i['id'][0] for i in hpo_helper.get_ancients([1]))
            #A = h[0] not in ancestor_cache[h[1]]
            #B = h[1] not in ancestor_cache[h[0]]
            #if A and B:
                #key = tuple(sorted([h[0],h[1]]))
                key = '-'.join(sorted([h[0],h[1]]))
                # using normaliser conveniently make the result consistently smaller than if h0 is a subclass of h1 if later using IC(h0)*IC(h1)/max_IC**2 for getting weights.
                normaliser = 2*min(len(hpo_freq[mode][h[0]]),len(hpo_freq[mode][h[1]]))
                result[key] += (method_cache[h[0]] + method_cache[h[1]]) / normaliser

    with open(phenopolis_utils.OFFLINE_CONFIG['hpo']['matrix_file'],'w') as outf:
        json.dump(result,outf)
