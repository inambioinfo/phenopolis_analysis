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
		Genon_sum: 64.4263097926
		Genon_hratio: 0.985077677075
		Genon_vratio: 0.992700643978
		Genon_sratio: 0.930102522342
		Genon_combined: 59.9230732433
	HP:0000707
		Genon_sum: 28.319483674
		Genon_hratio: 0.927743170062
		Genon_vratio: 0.896760844442
		Genon_sratio: 0.734355582504
		Genon_combined: 20.7965709296
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
		Genon_sum: 64.4263097926
		Genon_hratio: 0.985077677075
		Genon_vratio: 0.992700643978
		Genon_sratio: 0.930102522342
		Genon_combined: 59.9230732433
...
```
