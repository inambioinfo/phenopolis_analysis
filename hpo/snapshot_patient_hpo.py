#!/bin/env python
'''
this script is to snapshot patient hpo information for further analyses, such as gene-hpo relationship
'''
from __future__ import print_function, division
import sys
sys.path.append('../commons')
import phenopolis_utils
import os
import os.path

dbs = phenopolis_utils.get_mongo_collections()


unrelated_file = os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['unrelated_file'])
if os.path.isfile(unrelated_file):
    unrelated = open(unrelated_file,'r').readlines()
    unrelated = [i.rstrip() for i in unrelated]
else:
    unrelated = []

if '--out' in sys.argv:
    outfile=sys.stdout
    outf=outfile
else:
    outfile = os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['snapshot_file'])
    outf=open(outfile,'w')

outf.write('#p_id   unrelated   HPO\n')

for p in dbs['patient_db'].patients.find({},{'features':1, 'external_id':1}):
    #p_id   unrelated   hpo
    if 'features' not in p:
        continue
    unrelated_flag = 1 if p['external_id'] in unrelated else 0
    hpos = [i['id'] for i in p['features'] if i['observed'] == 'yes']
    if not hpos:
        continue
    # replace obsolete hpos
    hpos = [phenopolis_utils.replace_hpo(dbs['hpo_db'], [h,h])[0] for h in hpos]
    all_hpos = [] # union of all ancestors
    for hpo in hpos:
        anc = phenopolis_utils.get_hpo_ancestors(dbs['hpo_db'], hpo)
        for a in anc:
            all_hpos.extend(a['id'])
    all_hpos = set(all_hpos)
    #write to file
    outf.write('%s\t%s\t%s\n' % (p['external_id'],unrelated_flag,','.join(all_hpos)))


