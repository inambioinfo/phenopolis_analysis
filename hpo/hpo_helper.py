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
