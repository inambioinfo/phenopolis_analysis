from __future__ import print_function
import pymongo

conn = pymongo.MongoClient('phenotips.cs.ucl.ac.uk')
db=conn['uclex']
exac=conn['exac']

print('gene_name','gene_length','variants_num','pli',sep=',')
for gene in db.genes.find():
    gene_length=gene['stop']-gene['start']
    variants_num=len([v for v in db.variants.find({'canonical_gene_name_upper':gene['gene_name']})])
    if variants_num == 0: continue
    g=exac.pli.find_one({'gene':gene['gene_name']})
    if g: pli=g.get('pLI','NA')
    else: pli='NA'
    print(gene['gene_name'],gene_length,variants_num,pli,sep=',')
