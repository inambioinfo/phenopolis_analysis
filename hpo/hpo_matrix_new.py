#!/bin/env python
'''
generate weighted adjacency matrix for HPO [unrelated]
'''
from __future__ import print_function, division
import sys
sys.path.append('../commons')
import phenopolis_utils
import hpo_helper
import math
import json
import os
import pandas as pd
from collections import defaultdict
'''
main logic
'''
def build_matrices(snapshot_file,dbs,methods):
    print('read snapshot file')
    snapshot_df = pd.read_table(snapshot_file)
    print('remove patients with only HP:0000001 or HP:0000118')
    snapshot_df = snapshot_df[(snapshot_df['unrelated'] == 1) & ~snapshot_df['HPO'].isin(['HP:0000001','HP:0000118'])].reset_index()
    print('get unique term_set')
    term_set = get_term_set(snapshot_df,'HPO',',')
    print('build patient hpo df')
    df = build_patient_hpo_df(snapshot_df, 'HPO', term_set)
    # fill nan. not necessary if all patients have HPO term
    df = df.fillna(0)
    # get HPO freq
    df = df.astype(int)
    freq = df.sum()
    print('get IC')
    IC = IC_maker(freq['HP:0000001'])
    ic_df = freq.apply(IC)
    print('build co-oc matrix')
    df = df.T.dot(df)
    print('find all ancestors')
    ancestors = {}
    for i in freq.index:
        ancestors[i] = [j['id'][0] for j in phenopolis_utils.get_hpo_ancestors(dbs['hpo_db'],i)]
    print('build sim matrix')
    result_dfs = {}
    for k,v in methods.items():
        result_dfs[k] = v(dbs,df,ic_df,freq,ancestors)
    return result_dfs
'''
get unique terms
'''
def get_term_set(df,col,delimiter):
    def union(a,b):
        if pd.isnull(b):
            return a
        if not isinstance(b,str):
            return None
        return delimiter.join([a,b])
    term_set = df.apply(lambda x: reduce(union,x), axis=0)[col]
    return set(term_set.split(delimiter))

'''
build patient-hpo df
'''
def build_patient_hpo_df(df,col,term_set):
    result = pd.DataFrame()
    for t in term_set:
        result[t] =df[col].str.contains(t)
    return result

'''
make wrapper of hpo_helper.IC
'''
def IC_maker(t):
    def IC_wrapper(a):
        return hpo_helper.IC(a,t)
    return IC_wrapper

'''

'''
def beta(h_1, **kwargs):
    # in buffer?
    if kwargs['mode'] != 'asym':
        k = tuple(sorted([h_1['index'],kwargs['h2']]))
        if k in kwargs['buff']:
            return kwargs['buff'][k]
    ic1 = h_1[0]
    ic2 = kwargs['ic_df'][kwargs['h2']]
    freq1 = kwargs['freq'][h_1['index']]
    freq2 = kwargs['freq'][kwargs['h2']]
    if kwargs['mode'] == 'sym':
        result = 0.5 * (ic1 + ic2) / ( min(freq1,freq2) * kwargs['max_ic'] )
    elif kwargs['mode'] == 'asym': 
        result = 0.5 * (ic1 + ic2) / ( freq1 * kwargs['max_ic'] )
    elif kwargs['mode'] == 'lin':
        result = 2 * (kwargs['ic_df'][nca(h_1['index'],kwargs['h2'],kwargs['freq'],kwargs['ancestors'])]) / ( (ic1 + ic2) or 1 )
    elif kwargs['mode'] == 'resnik':
        result = kwargs['ic_df'][nca(h_1['index'],kwargs['h2'],kwargs['freq'],kwargs['ancestors'])]
    else:
        msg = 'Do not recognise similarity method: ' + kwargs['mode']
        raise ValueError(msg)
    if kwargs['mode'] != 'asym':
        kwargs['buff'][k] = result
    return result

'''
sym_WAM
'''
def sym_WAM(dbs,df,ic_df,freq,ancestors):
    print('build hpo_sym_WAM')
    weight_df = pd.DataFrame(index=ic_df.index)
    max_ic = hpo_helper.IC(1,freq['HP:0000001'])
    buff = {}
    for i,h in enumerate(ic_df.index):
        this = ic_df.reset_index().apply(beta,axis=1,dbs=dbs,ic_df=ic_df,freq=freq,max_ic=max_ic,h2=h,buff=buff,mode='sym')
        weight_df[h] = this.values
    return(df.multiply(weight_df))

'''
asym_WAM
'''
def asym_WAM(dbs,df,ic_df,freq,ancestors):
    print('build hpo_asym_WAM')
    weight_df = pd.DataFrame(index=ic_df.index)
    max_ic = hpo_helper.IC(1,freq['HP:0000001'])
    buff = {}
    for i,h in enumerate(ic_df.index):
        this = ic_df.reset_index().apply(beta,axis=1,dbs=dbs,ic_df=ic_df,freq=freq,max_ic=max_ic,h2=h,buff=buff,mode='asym')
        weight_df[h] = this.values
    return(df.multiply(weight_df))
   
'''
lin
'''
def lin(dbs,df,ic_df,freq,ancestors):
    print('build lin')
    result = pd.DataFrame(index=ic_df.index)
    buff = {}
    for i,h in enumerate(ic_df.index):
        print(i)
        this = ic_df.reset_index().apply(beta,axis=1,dbs=dbs,ic_df=ic_df,freq=freq,h2=h,buff=buff,mode='lin',ancestors=ancestors)
        result[h] = this.values
    return result

'''
resnik
'''
def resnik(dbs,df,ic_df,freq,ancestors):
    print('build resnik')
    result = pd.DataFrame(index=ic_df.index)
    buff = {}
    for i,h in enumerate(ic_df.index):
        print(i)
        this = ic_df.reset_index().apply(beta,axis=1,dbs=dbs,ic_df=ic_df,freq=freq,h2=h,buff=buff,mode='resnik',ancestors=ancestors)
        result[h] = this.values
    return result

'''
get nearest common ancestor
'''
def nca(h1,h2,freq,ancestors):
    ca = set(ancestors[h1]) & set(ancestors[h2])
    m = c = None
    for h in ca:
        if not m or freq[h] < c:
            c = freq[h]
            m = h
    return m
'''
main
'''
if __name__ == '__main__':
    snapshot_file = os.path.join('..',phenopolis_utils.OFFLINE_CONFIG['hpo']['snapshot_file'])
    dbs = phenopolis_utils.get_mongo_collections()
    methods = {
            'hpo_sym_WAM': sym_WAM,
            'hpo_asym_WAM': asym_WAM,
            'lin': lin,
            'resnik': resnik,
    }
    methods = {
            'lin': lin,
            'resnik': resnik,
    }
    store = pd.HDFStore('../data/public/hpo/store.h5')
    dfs = build_matrices(snapshot_file,dbs,methods)
    print('write to HDF5')
    for k,v in dfs.items():
        store[k] = v
    store.close()
    print('done')
