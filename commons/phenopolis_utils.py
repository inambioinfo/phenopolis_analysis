'''
some common functions used across different analyses
'''
from __future__ import print_function, division
import pymongo
import ConfigParser
import os
import errno
import sys
import sqlite3
import re
import itertools
import logging

'''
constants
'''
VALID_CHROMOSOMES = [str(i) for i in range(1,23)] + ['X','Y']

'''
parse config file, and make config global. If test, set DB_HOST as 'localhost'
'''
def _parse_config():
    # return {'section':{'key1':'value1'...},...}
    config = ConfigParser.ConfigParser()
    
    # get this path, and then read common.cfg
    path = os.path.dirname(os.path.abspath(__file__))
    config.read(os.path.join(path,'common.cfg'))
    result = {}
    for section in config.sections():
        options = config.options(section)
        result[section] = {}
        for option in options:
            result[section][option] = config.get(section, option)
    return result

OFFLINE_CONFIG = _parse_config()
# log to file
logging.basicConfig(filename=OFFLINE_CONFIG['debug']['log_file'],
        level=getattr(logging,OFFLINE_CONFIG['debug']['log_level'].upper()))
'''
get useful mongo collections
'''
def get_mongo_collections(test=None):
    if OFFLINE_CONFIG['mongodb']['db_port']:
        conn = pymongo.MongoClient(
            host = OFFLINE_CONFIG['mongodb']['db_host'],
            port = int(OFFLINE_CONFIG['mongodb']['db_port']),
        )
    else:
        conn = pymongo.MongoClient(
            host = OFFLINE_CONFIG['mongodb']['db_host'],
        )
    if not test:
        return {
            'hpo_db': conn[OFFLINE_CONFIG['mongodb']['db_name_hpo']],
            'phenopolis_db': conn[OFFLINE_CONFIG['mongodb']['db_name']],
            'patient_db': conn[OFFLINE_CONFIG['mongodb']['db_name_patients']],
            'pubmedbatch': conn[OFFLINE_CONFIG['mongodb']['db_name_pubmedbatch']],
        }

    else:
        return {
            'hpo_db': conn[test['hpo_db']],
            'phenopolis_db': conn[test['phenopolis_db']],
            'patient_db': conn[test['patient_db']],
            'pubmedbatch': conn[test['pubmedbatch_db']],
        }

'''
given chromosomes and db, return genes
'''
def get_chrom_genes(chroms, db):
    # give chrom numbers, get all genes on them
    result = []
    for chrom in chroms:
        chrom = str(chrom)
        if chrom not in VALID_CHROMOSOMES:
            raise ValueError('Error: %s is not a valid chromosome!' % chrom)
        genes = [g['gene_id'] for g in db.genes.find({'chrom':chrom})]
        result.extend(genes)
    return result

'''
given symbols, return ENSEMBL ids.
return as it is if one of the symbols already is ENSEMBL ids
Note that if a symbol cannot be found, it issues a warning.
'''
def symbols_to_ids(symbols,db):
    result = []
    ss = []
    fields = {
            '_id':0,
            'gene_id':1,
            'gene_name':1,
            }
    for s in symbols:
        if s.startswith('ENSG'):
            result.append(s)
        else:
            ss.append(s)
    this = db.genes.find({'gene_name':{'$in':ss}},fields)
    gene_ids = set(result + [i['gene_id'] for i in this])
    print(set([i['gene_name'] for i in this]))
    missing_symbols = set(ss) - set([i['gene_name'] for i in this])
    if missing_symbols:
        logging.warning('symbols not exist in the database when translating'+\
                ' to ENSEMBL ids: {}'.format(', '.join(missing_symbols)))
    return list(gene_ids)
'''
mkdir -p
http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
basically, it only makes one system call, therefore avoid racing problems.
not required for python3, can use `os.makedirs(name,mode,exist_ok=True)
'''
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

'''
some hpo ids are obsolete
this is a copy of lookups.py's replace_hpo. lookups.py is not working atm
'''
def replace_hpo(hpo_db, hpo):
    # some hpo_ids are obsolete.
    record = hpo_db.hpo.find_one({'id':hpo[0]})
    if not record:
        print ('no record in replace_hpo')
        print (hpo)
    if 'replaced_by' in record:
        new = hpo_db.hpo.find_one({'id':record['replaced_by'][0]})
        return [new['id'][0], new['name'][0]]
    else:
        return hpo

'''
get all ancestor nodes of a given hpo_id.
'''
def get_hpo_ancestors(hpo_db, hpo_id):
    """
    Get HPO terms higher up in the hierarchy.
    """
    h=hpo_db.hpo.find_one({'id':hpo_id})
    #print(hpo_id,h)
    if 'replaced_by' in h:
        # not primary id, replace with primary id and try again
        h = hpo_db.hpo.find_one({'id':h['replaced_by'][0]})
    hpo=[h]
    if 'is_a' not in h: return hpo
    for hpo_parent_id in h['is_a']:
        #p=hpo_db.hpo.find({'id':hpo_parent_id}):
        hpo+=list(itertools.chain(get_hpo_ancestors(hpo_db,hpo_parent_id))) 
    #remove duplicates
    hpo={h['id'][0]:h for h in hpo}.values()
    return hpo

'''
get common ancestors of two given hpos
'''
def get_hpo_common_ancestors(hpo_db, h1, h2):
    # return a list of hpo ids for h1 and h2's common ancestors
    a1 = get_hpo_ancestors(hpo_db, h1)
    a2 = get_hpo_ancestors(hpo_db,h2)
    an1 = []
    an2 = []
    for a in a1:
        an1.extend(a['id'])
    for a in a2:
        an2.extend(a['id'])
    return list(set(an1) & set(an2))

'''
minimise a list of hpos
'''
def hpo_minimum_set(hpo_db, hpo_ids=[]):
    '''
    minimize the hpo sets
    results = {'HP:0000505': [ancestors]}
    '''
    hpo_ids = list(set(hpo_ids))
    results = dict([(hpo_id, [ h['id'][0] for h in get_hpo_ancestors(hpo_db, hpo_id)],) for hpo_id in hpo_ids])
    # minimise
    bad_ids = []
    for i in range(len(hpo_ids)):
        for j in range(i+1,len(hpo_ids)):
            if hpo_ids[i] in results[hpo_ids[j]]:
                # i is j's ancestor, remove
                bad_ids.append(hpo_ids[i])
                break
            if hpo_ids[j] in results[hpo_ids[i]]:
                # j is i's ancestor, remove
                bad_ids.append(hpo_ids[j])
    return list(set(hpo_ids) - set(bad_ids))

'''
translate gene_names to ensembl ids. db = dbs['phenopolis_db']
'''
def gene_names_to_ids(db, queries):
    result = {}
    if not queries:
        return result
    gs = db.genes.find({'$or':
        [
            {'gene_name':{'$in':queries}},
            {'other_names':{'$in':queries}},
        ]
    })
    qs = set(queries)
    for g in gs:
        name = list(qs & set(g.get('other_names',[]) + [g['gene_name']]))[0]
        result[name] = {
                'id':g['gene_id'],
                'symbol':g['gene_name'],
                }

    return result


'''
get candidate genes and patients' hpos, solve, candidate genes, sex
'''
def get_candidate_genes(dbs, genes=None, fields=None):
    # set up some defaults. hpos = observed features.
    # solve would be 0 for unsolved and 1 for solved
    # sex 0 unknown, 1 male, 2 female
    # if genes == None, get all genes
    SEX_DICT = {
            'F': 2,
            'M': 1,
            'U': 0,
        }
    SOLVE_DICT = {
            'solved':1,
            'unsolved':0,
        }
    
    # fields of interests
    fields = fields or ['hpo','solve','genes','sex','external_id','contact']

    all_valid_p = [p for p in dbs['patient_db'].patients.find({}) if p.get('genes',[])]
    result = {}
    gene_names = []
    for k1 in all_valid_p:
        for k2 in k1['genes']:
            # there's one patient that has Somatic NLRP3 as gene.
            # and there's one patient has GPR98. should be ADGRV1
            if k2['gene'] == 'Somatic NLRP3':
                gene_names.append('NLRP3')
                continue
            if k2['gene'] == 'GPR98':
                gene_names.append('ADGRV1')
                continue
            # illegal char?
            k2['gene'] = k2['gene'].strip()
            if re.search(r'[^a-zA-Z0-9-]', k2['gene']):
                raise ValueError('Error: Illegal gene name "%s"' % k2['gene'])
            if not genes or k2['gene'] in genes:
                gene_names.append(k2['gene'])

    gene_dict = gene_names_to_ids(dbs['phenopolis_db'],gene_names)
    for p in all_valid_p:
        # deal with hpo and solve and sex
        temp  = {f:p.get(f,None) for f in fields}
        if 'hpo' in fields:
            temp['hpo'] = []
            for f in p['features']:
                if f['observed'] == 'yes':
                    f['id'],f['label'] = replace_hpo(dbs['hpo_db'],(f['id'],f['label']))
                    temp['hpo'].append(f)
        if 'solve' in fields:
            temp['solve'] = SOLVE_DICT[p['solved']['status']]
        if 'sex' in fields:
            temp['sex'] = SEX_DICT[p['sex']]
        for g in p['genes']:
            if not g['gene']: continue
            if g['gene'] == 'Somatic NLRP3':
                g['gene'] = 'NLRP3'
            if g['gene'] == 'GPR98':
                g['gene'] = 'ADGRV1'
            # when defined genes, gene_dict might not have g['gene']
            if g['gene'] not in gene_dict: continue
            gene_id = gene_dict[g['gene']]['id']
            result[gene_id] = result.get(gene_id,{
                'symbol':gene_dict[g['gene']]['symbol'],
                'data':[],
                })
            result[gene_id]['data'].append(temp)
    return result

