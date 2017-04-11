#!/usr/bin/env Rscript

# gene similarity based on the HPO terms
library(ontologyIndex)
library(ontologyPlot)
library(ontologySimilarity)
#set.seed(1)
hpo.ontology <- ontologyIndex::get_ontology('/cluster/scratch3/vyp-scratch2/reference_datasets/HPO/hp.obo')
#hpo.ontology <- ontologyIndex::get_ontology('wget http://purl.obolibrary.org/obo/hp.obo') 
#hpo.ontology <- ontologyIndex::get_ontology('hp.obo') 
d <- fread('/cluster/scratch3/vyp-scratch2/reference_datasets/HPO/ALL_SOURCES_TYPICAL_FEATURES_diseases_to_genes_to_phenotypes.txt',skip=1,data.table=FALSE)
colnames(d) <- c('diseaseId','symbol', 'geneid', 'hpo_id', 'hpo_name')

gene2hpo = list()
for (g in unique(d$symbol)) {
    gene2hpo[[g]] <- as.character(d[which(d$symbol==g),'hpo_id'])
}

x <- lapply(hpo, function(x_i) minimal_set(hpo.ontology, c("HP:0000001", x_i)))

sim_mat <- get_sim_grid(ontology=hpo.ontology, term_sets=x)

100*length(intersect(hpo2gene[[1]],hpo2gene[[2]]))/length(gene[[2]])

#map everyone's phenotype to non-redundant set
#information_content <- get_term_info_content(hpo, x)
#information_content <- descendants_IC(hpo)
#term_sets <- replicate(simplify=FALSE, n=7, expr=minimal_set(hpo, sample(hpo$id, size=8)))
#colnames(sim_mat) <- rownames(sim_mat) <- d$eid

range(as.numeric(sim_mat))

library(fastcluster)

h <- hclust(as.dist(sim_mat))


genes <- c('COX10', 'ECHS1', 'FOXRED1', 'LIPT1', 'MTFMT', 'ND3', 'NDUFA10', 'NDUFA12', 'NDUFA2', 'NDUFA4', 'NDUFA9', 'NDUFAF2', 'NDUFAF5', 'NDUFAF6', 'NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2', 'PDHA1', 'PET100', 'SLC19A3', 'SURF1', 'TACO1', 'TRNK', 'TRNV', 'TRNW')

pdf('~/hpo-cluster1.pdf')
onto_plot(hpo.ontology, term_sets=minimal_set(hpo.ontology,unlist(hpo[genes])))
dev.off()


#genes <- <list of hpo profiles>
#gene[c("ACTN1","ABCA4")]

#onto_plot(hpo, term_sets=genes[c("ACTN1","ABCA4")])
#cluster <- c("ACTN1","ABCA4")

get_sim_p(sim_mat, group=genes)

genes_with_ancestors <- lapply(hpo[genes], function(x) get_ancestors(hpo.ontology, x))

sort(table(unlist(genes_with_ancestors)))

sort(table(unlist(hpo[genes])))

library(proxy)


gene = list()
for (h in unique(d$hpo_name)) {
    gene[[h]] <- as.character(d[which(d$hpo_name==h),'symbol'])
}
#x <- simil(gene,method='Jaccard')
all.genes <- unique(unlist(gene))
X <- matrix(0,nrow=length(gene),ncol=length(all.genes))
colnames(X) <- all.genes
rownames(X) <- names(gene)
for (h in names(gene)) {
    X[h,gene[[h]]] <- 1
}



