#!/bin/env python
'''
connect to patients db, get hpo frequencies, write to hpo_freq.tsv
only consider unrelated individuals calculated by KING
'''
from __future__ import print_function, division
import json
import sys
sys.path.append('../commons')
import phenopolis_utils
import os
dbs = phenopolis_utils.get_mongo_collections()

outfile = os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['hpo_freq_file'])

infile = os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['snapshot_file'])

release = infile.split('.')[2].split('_')[-1]

with open(infile, 'r') as inf:
    result = {'related':{},'unrelated':{},'release':release}
    for row in inf:
        if row[0] == '#': continue
        row = row.rstrip().split('\t')
        p = row[0]
        unrelated = int(row[1])
        hpos = row[2].split(',')
        for h in hpos:
            # union
            result['related'][h] = result['related'].get(h,[])
            result['related'][h].append(p)
            if unrelated:
                result['unrelated'][h] = result['unrelated'].get(h,[])
                result['unrelated'][h].append(p)

    print('related patients length: %s' % len(result['related']['HP:0000001']))
    print('unrelated patients length: %s' % len(result['unrelated']['HP:0000001']))

with open(outfile,'w') as outf:
    json.dump(result,outf)
print('done')
