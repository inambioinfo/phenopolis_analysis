#!/usr/bin/env Rscript

library(magrittr)
library(ontologyIndex)
library(ontologySimilarity)
library(ontologyPlot)
#library(SimReg)

# human ontology
#hpo <- ontologyIndex::get_ontology('hp.obo',extract_tags='everything')
#old.hpo <- ontologyIndex::get_ontology('/cluster/scratch3/vyp-scratch2/reference_datasets/HPO/hp.obo')
#hpo.names <- unique(c(hpo$name,old.hpo$name))
#hpo.ids <- unique(c(hpo$ids,old.hpo$ids))

#setwd('/SAN/vyplab/NCMD/DECIPHER')
d <- read.csv('decipher-patients.csv')
decipher.phenotypes <- unlist(strsplit(d$phenotype,';'))
decipher.phenotypes <- decipher.phenotypes[which(decipher.phenotypes!='')]
X3 <- apply(d,1, function(x) {
    X2 <- rep(0,length(unique(decipher.phenotypes)))
    names(X2) <- unique(decipher.phenotypes)
    X2[intersect(unlist(strsplit(x[['phenotype']],';')),names(X2))] <- 1
    return(X2)
}) 
X3 <- t(X3) 
X <- as.matrix(X3)
out <- crossprod(X)  # Same as: t(X) %*% X
colnames(out) <- rownames(out) <- unique(decipher.phenotypes)


# prob of co-occurence
X4 <- colSums(X3)/nrow(X3)
X4 <- (X4%*%t(X4))
diag(X4) <- colSums(X3)/nrow(X3)
rownames(X4) <- colnames(X4)

head(colSums(X3))

head(sort(out[,'Nystagmus']/rowSums(out),decreasing=T),30)



exp.hpo.matrix <- matrix(0, nrow=length(unique(decipher.phenotypes)),ncol=length(unique(decipher.phenotypes)))


pheno.prob <- sort(table(decipher.phenotypes)/nrow(d),decreasing=TRUE)
colnames(exp.hpo.matrix) <- rownames(exp.hpo.matrix) <- names(pheno.prob)


for (n in names(pheno.prob)) {
exp.hpo.matrix[n,n]  <- pheno.prob[[n]]
}












setdiff(,hpo.names)


for (i in 1:nrow(d)) {

print(hpo.names[unlist(strsplit(d[i,'phenotype'],';'))])
break

}


