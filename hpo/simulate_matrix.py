#!/bin/env python
'''

'''
from __future__ import print_function, division
import sys
sys.path.append('../commons')
import phenopolis_utils
import os
import itertools
from optparse import OptionParser
from collections import defaultdict
import hpo_helper

'''
return the output file
'''
def get_outfile(input):
    outfolder = os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['simulation_matrix_folder'])
    # mkdir -p output folder
    phenopolis_utils.mkdir_p(outfolder)
    # get input file number
    basename = os.path.basename(input)
    num = basename.split('.')[0].split('_')[1]
    outfile = os.path.join(outfolder,'hpomatrix_'+num+'.json')
    return outfile

'''
generate matrix
'''
def matrix(data,freq):
    # get total and ICmax
    total = len(data)
    IC_max = hpo_helper.IC(1,total)
    result = defaultdict(float)
    method_cache = {}
    n = total
    for d in data:
        print(n)
        n -= 1
        hpos = d['hpo']

        # apparently each hpo co-occur with itself
        for h in hpos:
            key = '-'.join([h,h])

            method_cache[h] = method_cache.get(h,hpo_helper.IC(freq[h],total)/(IC_max))
            result[key] = result.get(key, method_cache[h])
        # check each combination of the hpos 
        for h in itertools.combinations(hpos,2):
            key = '-'.join(sorted([h[0],h[1]]))
            # using normaliser conveniently make the result consistently smaller than if h0 is a subclass of h1 if later using IC(h0)*IC(h1)/max_IC**2 for getting weights.
            normaliser = 2*min(freq[h[0]],hpo_freq[h[1]])
            result[key] += (method_cache[h[0]] + method_cache[h[1]]) / normaliser
    return result

'''
main
'''
if __name__ == "__main__":
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)

    parser.add_option("--input",
                      dest="input",
                      help="input sim file?")
    (options, args) = parser.parse_args()
    if not options.input:
        msg = 'has to specify input'
        raise ValueError(msg)
    # get dbs
    dbs = phenopolis_utils.get_mongo_collections()
    # get input data
    input_data = hpo_helper.get_json(options.input)
    # expand hpos
    hpo_helper.expand_hpo(input_data)
    # get hpo_freq
    hpo_freq = hpo_helper.get_hpo_freq(input_data)
    result = matrix(input_data, hpo_freq)
    # get output file
    outfile = get_outfile(options.input)
    hpo_helper.write_json(result, outfile)
    print('done')
