setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
library(gdata)

cBiomut <- read.xls(xls="/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/cBio_CNA_LGG_BRCA_GBM.xlsx", method="tab")

cBiomut <- read.delim(file="/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/cBio_CNA_LGG_BRCA_GBM.xlsx", sep="\t")

diseases <- c("BRCA", "GBM") #no CNA thresh for LGG

for(x in diseases){ 
  thresh <- read.delim(paste(".", x, "CNA", "all_thresholded.by_genes.txt", sep="/"), sep="\t")
  names <- grep("TCGA", colnames(thresh), fixed=T, value=T)
  newnames<- list()
    for(y in 1:length(names)){
      newnames[[y]] <- strsplit(names[y], split=".", fixed=T)
      newnames <- paste(newnames[[y]][[1]][c(1:3)],sep="", collapse="-")
    }
  colnames(thresh)[4:length(colnames(thresh))] <- unlist(newnames)
  
  matstart <- grep(x, colnames(cBiomut), fixed=T)
  cBio <- (cBiomut[,(matstart+1):(matstart+4)])
  colnames(cBio) <- c("Chr", "Cytoband", "Number", "Genes")
  cBio[,"Genes"] <- as.character(cBio[,"Genes"])
    for(i in 1:length(rownames(cBio))){
      if(!is.na(cBio[i,"Chr"])){next}
      cBio[(i-1), "Genes"] <- paste(cBio[i, "Genes"], cBio[(i-1), "Genes"], collapse=" ")
    }
  cBioFin <- cBio[!is.na(cBio[,"Chr"]),]
  
  index <- intersect(thresh[,"Cytoband"], cBioFin[,"Cytoband"])
  
  threshFin <- thresh[thresh[,"Cytoband"] %in% index,]
  save(threshFin, file=paste(".", x, "CNA/threshold_cyto_mat.RData", sep="/")) 
}


##########Trying to the the Genes column to count genes separately and not as groups
genes1 <- cBio[i, "Genes"]
genes2 <- cBio[i-1, "Genes"]
genesFin <- unlist(list(genes1, genes2))
blah <- as.character(paste(genesFin[1:length(genesFin)], collapse=","))
blah <- gsub(" ", ",", blah, fixed=T)


strsplit(as.character(genesFin)[,split=" ", fixed=T)
gsub(" ", "_", (as.character(genesFin)))

write.table(cBio, file="./Temp.txt", sep="\t")
temp <- read.delim("./Temp.txt", sep="\t")
scan(file="./Temp.txt")

genes1 <- temp[i, "Genes"]
genes2 <- temp[i-1, "Genes"]
genesFin <- unlist(list(genes1, genes2))
