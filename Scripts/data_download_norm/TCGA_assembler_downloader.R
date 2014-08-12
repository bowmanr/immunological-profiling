################ Loading TCGA ###############
setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler")

source("./Scripts/Module_A.r")
source("./Scripts/Module_B.r")
################ Clinical Data ################
master_folder <- "/Volumes/MacintoshHD3/data/TCGA_Assembler/Data/"
setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data/")
diseases <- list.files(recursive=F)

for(x in diseases){
#clinicalData <- DownloadClinicalData(traverseResultFile = "/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/DirectoryTraverseResult_Jan-30-2014.rda", saveFolderName = paste(master_folder, x, "/CD", sep = ""), cancerType = x, clinicalDataType = c("patient", "drug", "follow_up", "radiation"))

################ CNA Data ################
#CNAData <- DownloadCNAData(traverseResultFile = "/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/DirectoryTraverseResult_Jan-30-2014.rda", saveFolderName = paste(master_folder, x, "/CNA", sep = ""), cancerType = x, assayPlatform = "genome_wide_snp_6")

################ RNASeq Data ################
#RNASeq2Data <- DownloadRNASeqData(traverseResultFile = "/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/DirectoryTraverseResult_Jan-30-2014.rda", saveFolderName = paste(master_folder, x, "/RNASeq", sep = ""), cancerType = x, assayPlatform = "RNASeqV2", dataType = "rsem.genes.normalized_results")

############### Methylation Data #############
#MethylationData27 <- DownloadMethylationData(traverseResultFile = "/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/DirectoryTraverseResult_Jan-30-2014.rda", saveFolderName = paste(master_folder, x, "/Methylation", sep = ""), cancerType = x, assayPlatform = "humanmethylation27")
MethylationData450 <- DownloadMethylationData(traverseResultFile = "/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/DirectoryTraverseResult_Jan-30-2014.rda", saveFolderName = paste(master_folder, x, "/Methylation", sep = ""), cancerType = x, assayPlatform = "humanmethylation450")

############## RPPA Data ####################
#RPPAData <- DownloadRPPAData(traverseResultFile = "/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/DirectoryTraverseResult_Jan-30-2014.rda", saveFolderName = paste(master_folder, x, "/RPPA", sep = ""), cancerType = x, assayPlatform = "mda_rppa_core")
}