from __future__ import print_function
import pymongo

conn = pymongo.MongoClient('phenotips.cs.ucl.ac.uk')
db=conn['uclex']

print('gene_name','gene_length','variants_num',sep=',')
for gene in db.genes.find():
    gene_length=gene['stop']-gene['start']
    variants_num=len([v for v in db.variants.find({'canonical_gene_name_upper':gene['gene_name']})])
    if variants_num == 0: continue
    print(gene['gene_name'],gene_length,variants_num,sep=',')
