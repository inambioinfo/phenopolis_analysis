#!/bin/env python
from __future__ import print_function, division
import sys
import os
sys.path.append('../commons')
import phenopolis_utils

dbs = phenopolis_utils.get_mongo_collections()
outf = open(os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['generic']['patient_mini_file']), 'w')

with open(os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['generic']['patient_info_file']),'r') as inf:
    for row in inf:
        if row[0] == '#': 
            outf.write(row)
            continue
        row = row.rstrip().split('\t')
        hpos = phenopolis_utils.hpo_minimum_set(dbs['hpo_db'], hpo_ids=row[2].split(','))
        row[2] = ','.join(hpos)
        outf.write('\t'.join(row) + '\n')

