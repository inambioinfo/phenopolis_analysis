#!/bin/env python
from __future__ import print_function, division
import sys
import os
sys.path.append('../commons')
import phenopolis_utils
import math
import json
import copy
import random
from collections import defaultdict

'''
record the number of singletons.
they don't contribute anything to the final result anyway, like HP:0000001
'''
def get_singletons(data):
    # data is already reverse sorted on size, so don't worry
    result = defaultdict(int)
    good_hpos = set()
    rem_inds = []
    for d in data:
        if d['size'] > 1:
            good_hpos = good_hpos | set(d['hpo'])
        else:
            if d['hpo'][0] not in good_hpos:
                result[d['hpo'][0]] += 1
    return result

'''
wrapper of phenopolis ancestor
'''
def get_ancestors(db,id):
    return [i['id'][0] for i in phenopolis_utils.get_hpo_ancestors(db,id)]

'''
to see if a set is minimised
'''
def minimise_bool(dbs,hpos,new_hpo):
    for h in hpos:
        A = h in get_ancestors(dbs['hpo_db'], new_hpo)
        B = new_hpo in get_ancestors(dbs['hpo_db'],h)
        if A or B:
            return False
    return True

'''
get original snapshot. not interested in related
'''
def read_data():
    with open(os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['mini_snapshot_file']),'r') as inf:
        result = []
        for row in inf:
            if row[0] == '#': continue
            row = row.rstrip().split('\t')
            unrelated = int(row[1])
            if not unrelated: continue
            hpos = row[2].split(',')
            result.append( {'size':len(hpos), 'hpo':hpos} )
    return sorted(result, key=lambda x:x['size'], reverse=True)

def get_hpo_pool(data,singletons,singleton_cutoff):
    pool = []
    for d in data:
        if d['size'] > 1:
            pool.extend(d['hpo'])
        elif singletons[d['hpo'][0]] < singleton_cutoff:
            pool.extend(d['hpo'])
    return pool
    
'''
randomly generate given number of minimised set of hpos.
pool gets reduced and returned
'''
def get_random_hpos(dbs,pool,n):
    result = []
    if n == 1:
        chosen = random.sample(pool, 1)[0]
        result.append(chosen)
        pool.remove(chosen)
    else:
        for biblibaba in range(n):
            print(set(pool), result)
            while True:
                chosen = random.sample(pool, 1)[0]
                if minimise_bool(dbs, result, chosen):
                    result.append(chosen)
                    pool.remove(chosen)
                    break
    return result

'''
empty slot, remember slot counts, and randomly choose from hpos pool
has to follow rules:
    has to be minimised!
'''
if __name__ == '__main__':
    # set random seed
    random.seed(12345)
    dbs = phenopolis_utils.get_mongo_collections()
    singleton_cutoff = int(phenopolis_utils.OFFLINE_CONFIG['hpo']['singleton_cutoff'])
    # get original_data
    data = read_data()
    # remove singleton
    singletons = get_singletons(data)
    hpo_pool = get_hpo_pool(data,singletons,singleton_cutoff)
    s_count = phenopolis_utils.OFFLINE_CONFIG['hpo']['simulation_count']
    outfolder = phenopolis_utils.OFFLINE_CONFIG['hpo']['simulation_folder']
    # simulate
    for i in range(int(s_count)):
        print(i)
        print('===')
        outfile = os.path.join('..',outfolder,'snapsim_'+str(i)+'.json')
        # it's always smart to leave single hpo patients last, to avoid not being able to get minimised hpos. snapshot was already sorted.
        pool_copy = copy.copy(hpo_pool)
        for j in data:
            # singleton? do not change
            if j['size'] == 1 and singletons[j['hpo'][0]] >= singleton_cutoff:
                continue
            j['hpo'] = get_random_hpos(dbs,pool_copy,j['size'])
        
        with open(outfile,'w') as outf:
            json.dump(data,outf)
