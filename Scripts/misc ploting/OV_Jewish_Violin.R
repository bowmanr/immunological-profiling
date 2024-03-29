setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
library(matrixStats)
library(beeswarm)
library(caroline)
library(scales)
library(gdata)
###########################################
#Plotting LGG Gene Expression Based on Hist. Diagnosis
###########################################
folder <- "JewishHistory" #Name of folder for results
data <- "ssGSEA_Cyto_Module_RNASeqCounts.RData" #Name of data input

x <- "OV"
y <- tolower(x)
system(paste("mkdir", paste(".", x, "Plots", folder, sep="/"), sep=" "))

gsva_files<- grep(data, list.files(paste(".",x,"Results", "ssGSEA_Cyto", sep="/"),full.names=TRUE),value=TRUE)
load(gsva_files, verbose=T)

jewish <- read.delim(paste(".", x, "CD", paste("nationwidechildrens.org_clinical_patient_", y,".txt", sep=""), sep="/"))
jewish <- jewish[,c("bcr_patient_barcode", "jewish_religion_heritage_indicator")]
jewish <- jewish[-c(1:2),] #Removes 2 lines of heading
rownames(jewish) <- jewish[,1]

jewishFin <- as.matrix(jewish[,2])
jewishFin <- as.matrix(gsub("[Not Available]", "No Jewish Heritage", x=jewishFin[,1], fixed=T))
rownames(jewishFin) <- rownames(jewish)

index <- intersect(colnames(gsva_results), rownames(jewishFin))

variable_of_interest <- as.matrix(jewishFin[index,])
GSVA_results <- gsva_results[,index]

writeLines(paste("Plotting for", x, Sys.time()))

color_index <- round(seq.int(from=1,to=14,by=14/length(levels(factor(variable_of_interest)))))
cols <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "aquamarine4","#3288BD", muted("dodgerblue2"), "#5E4FA2", muted("blueviolet"))[color_index]

j <- 1

pdf(file=paste("./", x, "/Plots/", folder, "/AllJewish_violin.pdf",sep=""))
par(mfrow=c(2,2))

for(i in 1:length(rownames(GSVA_results))){
  violins(GSVA_results[i,], by=levels(factor(variable_of_interest[,j])), main=rownames(GSVA_results)[i], xlab="",ylab="Relative enrichment score", deciles=FALSE, connect="median", connectcol="black", names=rep(c(""),length(color_index)),drawRect=T, col=cols)
  axis(1, labels = FALSE)
  text(1:length(color_index), par("usr")[3], srt = 45, adj = 1, labels = levels(factor(variable_of_interest)), xpd=T, cex=0.9) 
  mtext(1, text = "Stage", line = 3.5)
}
dev.off()