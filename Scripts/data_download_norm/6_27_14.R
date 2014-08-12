################ Loading TCGA ###############
setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler")

source("./Scripts/Module_A.r")
source("./Scripts/Module_B.r")
################ CNA Data Processing ################
diseases <- list.dirs("./Data", recursive=F)
master_folder <- "/Volumes/MacintoshHD3/data/TCGA_Assembler/"

for(x in diseases){
CNAinput <- grep("nocnv_hg19", list.files(x, recursive=T, full.names=T), value = T)
RNAinput <- grep("rsem", list.files(x, recursive=T, full.names=T), value = T)

CNAData <- ProcessCNAData(inputFilePath =  CNAinput,
                           outputFileName = "processedCNA",
                          outputFileFolder = x,
                           refGenomeFile = "./Scripts/SupportingFiles/Hg19GenePosition.txt")

RNAData <- ProcessRNASeqData(inputFilePath =  RNAinput,
                          outputFileName = "processedRNA",
                          outputFileFolder = x,
                          dataType = "GeneExp",
                          verType = "RNASeqV2")
}