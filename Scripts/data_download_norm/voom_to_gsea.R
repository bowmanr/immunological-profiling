setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
library(phenoTest)
library(GSEABase)
library(GSVAdata) 
library(GSVA)
library(parallel)
library(edgeR)
###########################################################################
# Load in Gene Lists
###########################################################################
#Xue Modules
load("/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/Human Gene Modules/ModuleLists.RData", verbose=T)
#Xue Cytokine List
load("/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/CytoList_final_gene_sets.RData", verbose=T)
finalset <- unlist(res_set, recursive=F)
#Merged Cytokine and Mod List
mergedlist <- c(noBlank, finalset)

load("/Users/bowmanr/Projects/Imm_profile/Datasets/Gene_sets/Galon_Engler_list.RData", verbose=T)
mergedlist <- galon_engler
##################### Start Loop over Diseases ############################
diseases <- list.files(recursive=F)[-c(10,15)] #REMOVES STAD and Merged
for(x in diseases){
  system(paste("mkdir", paste(".", x, "Results", "GSEA", sep="/"), sep=" "))

  mut_files<-  grep("cBio_MutsigCV_BB",list.files(paste(".",x,"Mutation",sep="/"),full.names=TRUE),value=TRUE)
  if(length(mut_files)==2){
    mutation_1 <- read.delim(mut_files[1],sep="\t",skip=1,row.names=1,header=TRUE)
    mutation_1 <- mutation_1[,-dim(mutation_1)[2]]
    mutation_2 <- read.delim(mut_files[2],sep="\t",skip=1,row.names=1,header=TRUE)
    mutation_2 <- mutation_2[,-dim(mutation_2)[2]]
    mutation<-cbind(mutation_1,mutation_2)
    mutation <- mutation[,-dim(mutation)[2]]
  }
  
  if(length(mut_files)==1){
    mutation <- read.delim(mut_files,sep="\t",skip=1,row.names=1,header=TRUE)
    mutation <- mutation[,-dim(mutation)[2]]
  }
  
  if(length(mut_files)==0){
    print(paste("No Mutation data for",x,sep=" "))
    next
  }

  load(paste(".", x, "RNASeq", "RNASeqCounts.RData", sep="/"), verbose=T)


  index <- intersect(colnames(RNASeqFinal), rownames(mutation))
  mutation_set <- mutation[index,]
  RNASeqFinal <- RNASeqFinal[,index]
  gsea_results_final <- list()
  voomOutput <- voom(RNASeqFinal, design=NULL)
  for(i in 1:length(colnames(mutation_set))){ #loops over all mutations
    mut_interest<-mutation_set[,i] 
    mut_interest <- ifelse(mut_interest=="NaN","Wild-Type","Mutant")
    if(length(levels(factor(mut_interest)))!=2){
        gsea_results_final[[i]] <- "NULL"
      next}
      design <- model.matrix(~0+factor(mut_interest))
      colnames(design) <- c("Mutant","Wild_type")
      
      mergedlist_final <- lapply(mergedlist,intersect,rownames(voomOutput$E))

      gsea_results_final[[i]] <- camera(y=voomOutput$E,index=mergedlist_final,design=design)

    }
  names(gsea_results_final) <- colnames(mutation_set)
  save(gsea_results_final,file=paste(".", x, "Results/GSEA", "immune_GSEA_Results.RData", sep="/")) 
}
