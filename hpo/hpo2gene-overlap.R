#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(optparse, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))

library(ontologyIndex)

option_list <- list( 
make_option(c("--hpo"), help = "hpo term")
)
OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

hpo.ontology <- ontologyIndex::get_ontology('/cluster/scratch3/vyp-scratch2/reference_datasets/HPO/hp.obo')
print(hpo.id <- opt$hpo)

hpo.mapping <- read.csv('/SAN/vyplab/NCMD/NikolasPontikos/hpo-id-name.csv')
#print(hpo.term <- hpo.mapping[which(hpo.mapping$hpo_id==hpo.id),'hpo_name'])

hpo2gene <- readRDS('/SAN/vyplab/NCMD/NikolasPontikos/hpo2gene.rds')


print(X <- hpo2gene[[hpo.id]])
print(length(X))
print(length(unique(X)))

d <- lapply(hpo2gene, function (x) {
    length(intersect(X,x))/length(X)
})

as.numeric(d) -> d

saveRDS(d,file=sprintf("dist/%s.RDS",hpo.id))



