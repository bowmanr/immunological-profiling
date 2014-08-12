setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
library(matrixStats)
library(beeswarm)
library(caroline)
library(scales)
library(gdata)
###########################################
#Plotting LGG Gene Expression Based on Hist. Diagnosis
###########################################
folder <- "MenopauseStatus" #Name of folder for results
data <- "ssGSEA_Cyto_Module_RNASeqCounts.RData" #Name of data input

x <- "BRCA"
y <- tolower(x)
system(paste("mkdir", paste(".", x, "Plots", folder, sep="/"), sep=" "))

gsva_files<- grep(data, list.files(paste(".",x,"Results", "ssGSEA_Cyto", sep="/"),full.names=TRUE),value=TRUE)
load(gsva_files, verbose=T)

meno <- read.delim(paste(".", x, "CD", paste("nationwidechildrens.org_clinical_patient_", y,".txt", sep=""), sep="/"))
meno <- meno[,c("bcr_patient_barcode", "menopause_status")]
meno <- meno[-c(1:2),] #Removes 2 lines of heading
rownames(meno) <- meno[,1]
menoFin <- as.matrix(meno[,2])
colnames(menoFin) <- "Status"
rownames(menoFin) <- rownames(meno)

noNA <- grep("[Not Available]", menoFin[,"Status"], value=T, fixed=T, invert=T)
noNA <- grep("[Not Evaluated]", noNA, value=T, fixed=T, invert=T)
noNA <- grep("[Unknown]", noNA, value=T, fixed=T, invert=T)
noNA <- grep("Indeterminate", noNA, value=T, fixed=T, invert=T)
noindex <- menoFin[,"Status"]==levels(factor(noNA))

menoFin <- as.matrix(menoFin[,"Status"][noindex])
menoFin <- as.matrix(gsub("Peri (6-12 months since last menstrual period)", "Peri", x=menoFin[,1], fixed=T))
menoFin <- as.matrix(gsub("Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)", "Post", x=menoFin[,1], fixed=T))
menoFin <- as.matrix(gsub("Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)", "Pre", x=menoFin[,1], fixed=T))

index <- intersect(colnames(gsva_results), rownames(menoFin))


variable_of_interest <- as.matrix(menoFin[index,])
GSVA_results <- gsva_results[,index]

writeLines(paste("Plotting for", x, Sys.time()))

color_index <- round(seq.int(from=1,to=14,by=14/length(levels(factor(variable_of_interest)))))
cols <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "aquamarine4","#3288BD", muted("dodgerblue2"), "#5E4FA2", muted("blueviolet"))[color_index]

j <- 1

pdf(file=paste("./", x, "/Plots/", folder, "/AllHist_violin.pdf",sep=""))
par(mfrow=c(2,2))

for(i in 1:length(rownames(GSVA_results))){
  violins(GSVA_results[i,], by=levels(factor(variable_of_interest[,j])), main=rownames(GSVA_results)[i], xlab="",ylab="Relative enrichment score", deciles=FALSE, connect="median", connectcol="black", names=rep(c(""),length(color_index)),drawRect=T, col=cols)
  axis(1, labels = FALSE)
  text(1:length(color_index), par("usr")[3]-0.03, srt = 45, adj = 1, labels = levels(factor(variable_of_interest)), xpd=T, cex=0.9) 
  mtext(1, text = "Stage", line = 3.5)
}
dev.off()