library(edgeR)
library(Rsubread)
library(matrixStats)
library(data.table)
library(reshape)

setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data/")
diseases <- list.files(recursive=F)[-14] #No RNASeq Data for STAD

for(x in diseases){
load(paste(".", x, "RNASeq/processed_RNA.rda", sep = "/"), verbose=T)

RNASeq <- Data
rownames(RNASeq) <- Des[,1]
genes <- Des[,1]!="?"
noQ <- Des[genes,1]

if(dim(as.matrix(table(table(noQ))))[1]==2){

dup_genes <- duplicated(noQ) | duplicated(noQ, fromLast = TRUE)
aggNames <- noQ[dup_genes]
aggInput <- (RNASeq)[aggNames,]

writeLines(paste("Aggregating genes for", x, Sys.time()))

RNASeqAgg <- aggregate(aggInput, by = list(aggNames), FUN = "mean")

rownames(RNASeqAgg) <- RNASeqAgg[,"Group.1"]
aggRNA <- RNASeqAgg[,-1]
RNAnoQ <- RNASeq[noQ,]
noQnoDup <- RNAnoQ[!dup_genes,]

idBreakdown <- do.call(rbind,strsplit(colnames(noQnoDup), split="-"))[,1:3]
newID <- apply(idBreakdown, 1, paste, collapse="-")
duplicated_IDS <- table(newID)
duplicated_index<-as.character(newID) %in% names(duplicated_IDS)[duplicated_IDS > 1]
singlets <- noQnoDup[,!duplicated_index]

writeLines(paste("Aggregating samples for", x, Sys.time()))
doublets <-noQnoDup[,duplicated_index] ; doublet_Ids <- newID[duplicated_index]

agg_doublets<-aggregate(t(noQnoDup[,duplicated_index]),by=list(newID[duplicated_index]),FUN="mean")
rownames(agg_doublets) <- agg_doublets[,1]
agg_doublets <- t(agg_doublets[,-1])
RNASeqFinal<- cbind(singlets,agg_doublets)

colnames(RNASeqFinal)<- apply(do.call(rbind,strsplit(colnames(RNASeqFinal), split="-"))[,1:3],1,function(x){paste(x,collapse="-")})

if(any(duplicated(colnames(RNASeqFinal)))) {
  print(paste("There are duplicates in the colnames of",x,sep=" "))
  next
}

save(RNASeqFinal, file=paste(".", x, "RNASeq", "RNASeqCounts.RData", sep="/"))
voomOutput <- voom(RNASeqFinal, design=NULL)
exprs_log2 <- voomOutput$E

exprs_zScore <- (exprs_log2-(rowMeans(exprs_log2)))/rowSds(exprs_log2)

rm(Data, Des)

save(exprs_zScore, file = paste(".", x, "RNASeq", "zScore.RData", sep="/"))
writeLines(paste("Saved all files for", x, Sys.time()))

}else{
  writeLines(paste(x,"has a weird duplicate things we need to look at
                    We need to write a special exception"))
}
}