#!/bin/env python
'''
get real hpo_freq data, and minify it. only interested in unrelated
'''
from __future__ import print_function, division
import json
import sys
sys.path.append('../commons')
import phenopolis_utils
import hpo_helper
import os

'''
minify
'''
def minify(data):
    result = {}
    for k,v in data['unrelated'].items():
        result[k] = len(v)
    return result

if __name__ == '__main__':
    # get infile and outfile
    infile = os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['hpo_freq_file'])
    outfile = os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['mini_real_freq_file'])
    # read infile
    data = hpo_helper.get_json(infile)
    # minify
    result = minify(data)
    # write outfile
    hpo_helper.write_json(result,outfile)
