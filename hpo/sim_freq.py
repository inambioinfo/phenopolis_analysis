#!/bin/env python
'''
test hpo frequency, and compare with the real one
'''
from __future__ import print_function, division
import sys
sys.path.append('../commons')
import phenopolis_utils
import os
import itertools
from optparse import OptionParser
import hpo_helper

'''
return the output file
'''
def get_outfile(input):
    outfolder = os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['simulation_freq_folder'])
    # mkdir -p output folder
    phenopolis_utils.mkdir_p(outfolder)
    # get input file number
    basename = os.path.basename(input)
    num = basename.split('.')[0].split('_')[1]
    outfile = os.path.join(outfolder,'hpofreq_'+num+'.json')
    return outfile

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
    # get dbs
    dbs = phenopolis_utils.get_mongo_collections()
    # get input data
    input_data = hpo_helper.get_json(options.input)
    # get total and ICmax
    total = len(input_data)
    IC_max = hpo_helper.IC(1,total)
    # expand hpos
    hpo_helper.expand_hpo(input_data)
    # get hpo_freq
    hpo_freq = hpo_helper.get_hpo_freq(input_data)
    # get output file
    outfile = get_outfile(options.input)
    hpo_helper.write_json(hpo_freq, outfile)
    print('done')
