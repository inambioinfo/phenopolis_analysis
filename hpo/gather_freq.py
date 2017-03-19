#!/bin/env python
'''
group all freq.json into one
'''
from __future__ import print_function, division
import sys
sys.path.append('../commons')
import phenopolis_utils
import os
import hpo_helper
from collections import defaultdict

'''
get input folder
'''
def get_folder():
    folder = phenopolis_utils.OFFLINE_CONFIG['hpo']['simulation_folder']
    folder = os.path.join('..',folder,'freq')
    return folder

'''
group data
'''
def group(folder):
    result = defaultdict(list)
    for f in os.listdir(folder):
        f = os.path.join(folder, f)
        data = hpo_helper.get_json(f)
        for k,v in data.items():
            result[k].append(v)
    return result

'''
main
'''
if __name__ == "__main__":
    # get input folder
    folder = get_folder()
    # group data
    group_data = group(folder)
    outfile = os.path.join('..', phenopolis_utils.OFFLINE_CONFIG['hpo']['simulation_freq_file'])
    hpo_helper.write_json(group_data, outfile)
    print('done')
