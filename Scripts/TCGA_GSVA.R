setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
library(GSEABase)
library(GSVAdata) 
library(GSVA)
library(parallel)
library(marray)
library(ggplot2)
################## gset.idx Possibilities ##########################
#Xue Modules
load("/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/Human Gene Modules/ModuleLists.RData", verbose=T)
#Galon-Engler
load("/Users/bowmanr/Projects/Imm_profile/Datasets/Gene_sets/Galon_Engler_list.RData", verbose=T)
#Xue Cytokine List
load("/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/misc data/CytoList_final_gene_sets.RData", verbose=T)
finalset <- unlist(res_set, recursive=F)
#Merged Cytokine and Mod List
mergedlist <- c(noBlank, finalset)
################## All Disesase GSVA/ssGSEA ##########################
modlist <- mergedlist #See above for options
method <- "ssgsea" #Method in GSVA
folder <- "ssGSEA_Cyto" #New folder made for results
name <- "ssGSEA_Cyto_Module_RNASeqCounts.RData" #Name of results, disease added in save function
diseases <- list.files(recursive=F)[-c(10,15)] #REMOVES STAD and Merged

for(x in diseases){
  load(paste(".", x, "RNASeq", "RNASeqCounts.RData", sep="/"), verbose=T)
  writeLines(paste("Loaded", x, sep=" "))
  
  writeLines(paste("Working on GSVA for", x, Sys.time()))
  system(paste("mkdir",paste(".",x,"Results", folder, sep="/"),sep=" "))
  
  gsva_results <- gsva(expr=as.matrix(RNASeqFinal), gset.idx.list=modlist, method=method, rnaseq=TRUE,
                       min.sz=0, max.sz=10000, parallel.sz = 10, verbose=TRUE)
  save(gsva_results, file = paste(".", x,"Results", folder, paste(x, name, sep=""), sep = "/"))
}
#Script runs into memory error if run in RStudio, no error in Terminal. Data was the same between the two, plotting commenced

################GSVA on Normal, non-TCGA datasets ##########################
setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Diseases_non_TCGA")

modlist <- mergedlist #See above for options
method <- "ssgsea" #Method in GSVA

diseases <- c("Brain", "Skin", "Breast", "Lung")

for(x in diseases){
  GSE <- grep("GSE", list.files(paste(".", x, sep="/"), full.names=T), value=T)
  for(y in 1:length(GSE)){
        exprs <- grep("exprs", list.files(GSE[y], full.names=T), value=T)
        if(length(exprs)==0){next}
        

        load(exprs, verbose=T);
        if(exprs=="GSE14905"){
          exprs_1 <- exprs_2
        }
        exprs_1 <- exprs_1[!apply(exprs_1,1,function(x) any(is.na(x))),]
        name <- strsplit(GSE[y], split="/")
        name <- paste(name[[1]][3], ".RData", sep="")
        
        gsva_results <-gsva(expr=as.matrix(exprs_1), gset.idx.list=modlist, method=method, rnaseq=TRUE, min.sz=0, max.sz=10000, parallel.sz = 10, verbose=TRUE)
        save(gsva_results, file = paste(".", "GSVA_Results", paste(x, name, sep="_"), sep="/"))
    }
}
################GSVA on merged datasets ##########################
setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data/Merged")

merged_files <- grep("RNASeqCounts",list.files(,recursive=F),value=TRUE)
for(x in merged_files){
  load(x, verbose=T)
  system(paste("mkdir",paste(".","Results", folder, sep="/"),sep=" "))
  
  writeLines(paste("Working on ssGSEA for", x, Sys.time()))
  gsva_results <- gsva(expr=as.matrix(bigmat), gset.idx.list=modlist, method=method, rnaseq=TRUE,
                       min.sz=0, max.sz=10000, parallel.sz = 24, verbose=TRUE)
  save(gsva_results, file = paste("./Results", folder, paste(x, names, sep = "_"), sep="/"))
}


in BRCA, COAD, GBM, HNSC, KIRC, KIRP, LGG, LUAD, LUSC, OV, READ, SARC, SKCM, THCA: rsession(66533) malloc: *** error for object 0x7ffa9aa01e80: pointer being freed was not allocated
*** set a breakpoint in malloc_error_break to debug



################Pheno for merged datasets ##########################

setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")

diseases <- list.files(recursive=F)[-c(10,15)] #REMOVES STAD
samples <-list()
for(x in 1:length(diseases)){
  load(paste(".", diseases[x], "RNASeq", "zScore.RData", sep="/"), verbose=T)
  writeLines(paste("Loaded", diseases[x]))
  samples[[x]] <- colnames(exprs_zScore)
}
names(samples) <- diseases

pheno<-lapply(samples,function(x){cbind(x,rep("disease",length(x)))}

pheno <- list()
for(i in 1:length(samples)){
	pheno[[i]]<-cbind(samples[[i]],rep(names(samples)[i],length(i)))
}
All_samples<-do.call(rbind,pheno)
save(All_samples,file="All_samples_pheno.RData")