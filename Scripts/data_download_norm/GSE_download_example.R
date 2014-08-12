source("/Users/bowmanr/Projects/generalScripts/affy_utils_1.R")
setwd("/Users/bowmanr/Projects/Imm_profile/Datasets/Diseases_non_TCGA")

pheno <- read.delim("GSE2361_series_matrix.txt",skip=45,nrows=76-29,header=TRUE)
data <- read.delim("GSE2361_series_matrix.txt",skip=76,header=TRUE,row.names=1)
exprs <- log2(data)

exprs_1<-summarize.exprs.data(exprs,rownames="hgu133a.db") #either if GPL96 hgu133a.db or GPL570 hgu133plus2.db

pheno_1 <- t(pheno[c(1,11),])
pheno_1[,2] <- do.call(rbind,strsplit(pheno_1[,2],split="Normal human tissue:Normal "))[,2]
rownames(pheno_1) <- pheno_1[,1]

save(exprs_1, file="GSE2361_exprs.RData")
save(pheno_1, file="GSE2361_pheno.RData")