################ Loading TCGA ###############
setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler")

source("./Scripts/Module_A.r")
source("./Scripts/Module_B.r")
################ CNA Data Processing ################
master_folder <- "/Volumes/MacintoshHD3/data/TCGA_Assembler/Data/"
setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data/")
diseases <- list.files(recursive=F)[-1]

for(x in diseases){
  writeLines("**********************************************************************************");
  writeLines("");
  writeLines(paste("Working on ", x,  sep = ""));
  
CNAinput <- grep("nocnv_hg19", list.files(paste(".",x,"CNA",sep="/"), full.names=T), value = T)
RNAinput <- grep("illuminahiseq_rnaseqv2__rsem.genes", list.files(paste(".",x,"RNASeq",sep="/"), full.names=T), value = T)

CNAData <- ProcessCNAData(inputFilePath =  CNAinput,
                          outputFileName = "processed_CNA",
                          outputFileFolder = paste(".",x,"CNA",sep="/"),
                          refGenomeFile = "../Scripts/SupportingFiles/Hg19GenePosition.txt")

RNAData <- ProcessRNASeqData(inputFilePath =  RNAinput,
                         outputFileName = "processed_RNA",
                          outputFileFolder = paste(".",x,"RNASeq",sep="/"),
                          dataType = "GeneExp",
                          verType = "RNASeqV2")
}