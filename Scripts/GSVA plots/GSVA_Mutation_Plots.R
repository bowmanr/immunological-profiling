setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
library(matrixStats)
library(beeswarm)
library(caroline)
library(scales)
library(gdata)
###########################################
#Plotting all Gene Expression Based on Mutation
###########################################
folder <- "Significant_Mutation_ssGSEA" #Name of folder for results
data <- "ssGSEA_Cyto_Module_RNASeqCounts.RData" #Name of data input minus the disease
cut <- 0.1 #pvalue cut-off for significance

diseases <- list.files(recursive=F)[-c(10,15)] #REMOVES STAD and Merged
for(x in diseases){
  system(paste("mkdir", paste(".", x, "Plots", folder, sep="/"), sep=" "))
  
  mut_files<-  grep("cBio_MutsigCV_BB",list.files(paste(".",x,"Mutation",sep="/"),full.names=TRUE),value=TRUE)
  
  gsva_files<- grep(data, list.files(paste(".",x,"Results", "ssGSEA_Cyto", sep="/"),full.names=TRUE),value=TRUE)
    ##Import and merge mutation files
  if(length(mut_files)==0){
    #print(paste(x,"does not have mutation data",sep=" "))
    next}
  
  if(length(gsva_files)==0){
    #print(paste(x,"does not have no GSVA Results",sep=" "))
    next}
  
  if(length(mut_files)==2){
    #print(paste(x,"has two mutation files",sep=" "))
    mutation_1 <- read.delim(mut_files[1],sep="\t",skip=1,row.names=1,header=TRUE)
    mutation_1 <- mutation_1[,-dim(mutation_1)[2]]
    mutation_2 <- read.delim(mut_files[2],sep="\t",skip=1,row.names=1,header=TRUE)
    mutation_2 <- mutation_2[,-dim(mutation_2)[2]]
    mutation<-cbind(mutation_1,mutation_2)
    mutation <- mutation[,-dim(mutation)[2]]
  }
  
  if(length(mut_files)==1){
    #print(paste(x,"has one mutation file",sep=" "))
    mutation <- read.delim(mut_files,sep="\t",skip=1,row.names=1,header=TRUE)
    mutation <- mutation[,-dim(mutation)[2]]
  }
  
  ##Import and merge mutation and GSVA results files
  load(gsva_files, verbose=T)
  
  index <- intersect(colnames(gsva_results), rownames(mutation))
  
  mutation_set <- mutation[index,]
  GSVA_results <- gsva_results[,index]
  
  writeLines(paste("Plotting for", x, Sys.time()))
  
  for(j in 1:length(colnames(mutation_set))){ #loops over all mutations
    mut_interest<-mutation_set[,j] 
    mut_interest <- ifelse(mut_interest=="NaN","Wild-Type","Mutant")
    
    if(length(levels(factor(mut_interest)))==1){
      #print(paste("No interesting mutation data for",names(variable_of_interest)[j],"in",x,sep=" "))
      next
    }
    if(length(levels(factor(mut_interest)))!=1){
      #print(paste("Lots of interesting mutation data for",names(variable_of_interest)[j],"in",x,sep=" "))
      
      if(any(as.matrix(table(mut_interest))[,1]<3)){
        next
      }
      
      pvalue <- p.adjust(apply(GSVA_results,1,function(x)t.test(x~factor(mut_interest))$p.value), n=length(rownames(GSVA_results)), method="fdr")
      
      if(any(pvalue<=(cut))){
        
        pdf(file=paste("./", x, "/Plots/", folder, "/", colnames(mutation_set)[j],".pdf",sep=""))
        par(mfrow=c(2,2))
        for(i in 1:length(rownames(GSVA_results))){ #loops over modules and cytokines
            if(pvalue[i]<=(cut)){
              
              boxplot(GSVA_results[i,]~factor(mut_interest),
                      main=rownames(GSVA_results)[i],
                      xlab="Variable of interest",ylab="Relative enrichment score",
                      names=levels(factor(mut_interest)), cex.axis=0.9,las=0)
              beeswarm(GSVA_results[i,]~factor(mut_interest),
                       method="hex", col=alpha("black", 0.5), corral=c("omit"), pch=19, add=TRUE)
            }
        }
        dev.off()
      }
    }
  }
}