library(edgeR)
library(Rsubread)
library(matrixStats)
library(data.table)
library(reshape)

setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data/")
diseases <- list.files(recursive=F)[-c(10, 15)] #Removes STAD and Merged
system(paste("mkdir", paste(".", "Merged", sep="/"), sep=" "))

Data <- list()

for(x in diseases){
  load(paste(".", x, "RNASeq", "RNASeqCounts.RData", sep = "/"), verbose=T)
  Data[[x]] <- RNASeqFinal
}
bigmat <- do.call(cbind,list(Data$COAD, Data$READ))# Edited manually to combine different diseases: LUSC&LUAD, KIRP&KIRC, GBM&LGG, COAD&READ
  
save(bigmat, file=paste(".", "Merged", "COAD_READRNASeqCounts.RData", sep="/"))
    voomOutput <- voom(bigmat, design=NULL)
    exprs_log2 <- voomOutput$E
    
    exprs_zScore <- (exprs_log2-(rowMeans(exprs_log2)))/rowSds(exprs_log2)
        
    save(exprs_zScore, file = paste(".", "Merged", "COAD_READzScore.RData", sep="/"))
