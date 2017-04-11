'''
convert decipher-patients data into hpo co-occurrence matrix
'''
from __future__ import print_function, division
import pandas as pd
from optparse import OptionParser
import sys

'''
convert csv to cooc matrix, main logic
'''
def csv_to_cooc(csv_file):
    # read csv
    print('read csv')
    csv_df = pd.read_csv(csv_file)
    # get unique terms
    print('get term set')
    term_set = get_term_set(csv_df)
    # build patient-hpo matrix
    print('build patient-hpo matrix')
    df = build_patient_hpo_df(csv_df, term_set)
    # fill nan
    df = df.fillna(0)
    # co-oc matrix
    print('build co-oc matrix')
    df = df.astype(int)
    return df.T.dot(df)

'''
get unique terms
'''
def get_term_set(df):
    def union(a,b):
        if pd.isnull(b):
            return a
        if not isinstance(b,str):
            return None
        return ';'.join([a,b])
    term_set = df.apply(lambda x: reduce(union,x), axis=0)['phenotype']
    return set(term_set.split(';'))

'''
build patient-hpo df
'''
def build_patient_hpo_df(df, term_set):
    result = pd.DataFrame()
    for t in term_set:
        result[t] =df['phenotype'].str.contains(t)
    return result


if __name__ == '__main__':
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)

    parser.add_option("--input",
                      dest="input", default='downloads/decipher-patients-sample.csv',
                      help="config everything there [default: %default]")
    (options, args) = parser.parse_args()

    result = csv_to_cooc(options.input)
    store = pd.HDFStore('output/store.h5')
    store['hpo_cooc'] = result
    print('done')
