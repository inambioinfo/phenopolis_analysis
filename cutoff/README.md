__Use Phenogenon to assess gene-phenotype relationships__

## Examples

1. When you have a set of phenotypes in mind, and know what mode of inheritance you are looking at:
```
from __future__ import print_function, division
from Genon import *
G = Genon()
G.init_genes()
GR = GenonResult()
perfect_result = G.analyse(GR)
print(perfect_result)
```
You will get (test data is needed)
```
SCN1A - Given mode is d
      - Number of patients with rare variants: 100
      - HPOs are given...
	HP:0001250
		HGF: 64.4263097926
		cadd_15_ratio: 0.992700643978
		MOI_score: 59.9230732433
	HP:0000707
		HGF: 28.319483674
		cadd_15_ratio: 0.896760844442
		MOI_score: 20.7965709296
...
```
2. When you just want to know: given a gene, what phenotypes and mode of inheritance are assocated with it?
```
G = Genon()
GR = GenonResult()
G.init_genes()
for g in G.genes:
    G.genes[g].hpos = None
    G.genes[g].mode = None
predict_both = G.analyse(GR)
print(predict_both)
```
You will get
```
SCN1A - Predicted mode is d
      - Number of patients for mode d: 100
      - Number of patients for mode r: 13
      - HPOs are predicted...
	HP:0001250
		HGF: 64.4263097926
		cadd_15_ratio: 0.992700643978
		MOI_score: 59.9230732433
...
```
