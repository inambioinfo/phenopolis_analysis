


```python
 python gene_analysis.py > gene_stats.txt
```

Open R:
```R
read.csv('gene_stats.txt')->d
plot(d$gene_length,d$variants_num,col=ifelse(d$pli>.9,'red','black'),pch=20)
```

