
setwd('/SAN/vyplab/NCMD/NikolasPontikos/dist/')

#hpo.dist<-as.matrix(do.call('rbind',lapply(list.files(pattern='*.RDS'),function(x){readRDS(x)})))
hpo.dist<-as.matrix(do.call('rbind',lapply(n$hpo_id,function(x){readRDS(paste(x,'RDS',sep='.'))})))

rownames(hpo.dist) <- colnames(hpo.dist)


