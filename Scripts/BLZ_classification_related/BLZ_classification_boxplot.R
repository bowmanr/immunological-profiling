setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
library(matrixStats)
library(beeswarm)
library(scales)
library(gdata)
###########################################
#Plotting BLZ/Vehicle
###########################################
diseases <- list.files(recursive=F)[-c(10,15)] #Removes STAD and Merged

for(x in diseases){

load(paste(".", x, "Results", "LassoRegression_Mouse_BLZ945.RData", sep="/"), verbose=T)
class <- as.matrix(exprs.predict[,"Class"])

load(paste(".", x, "Results", "RNASeq_counts_ssgsea_macrophage_modules.RData", sep="/"), verbose=T)
results <- gsva_results

index <- intersect(colnames(results), rownames(class))
print(paste(x,length(index)))

variable_of_interest <- as.matrix(class[index,])
GSVA_results <- results[,index]

writeLines(paste("Plotting for", x, Sys.time()))

pdf(file=paste("./", x, "/Plots/", "ssGSEA_BLZClass",".pdf",sep=""))
par(mfrow=c(2,2))
for(i in 1:length(rownames(GSVA_results))){ #this loops over each module
  boxplot(GSVA_results[i,]~factor(variable_of_interest), main=rownames(GSVA_results)[i],
          xlab="Variable of interest", ylab="Relative enrichment score",
          names=levels(factor(variable_of_interest)), cex.axis=0.9,las=0, cex.names=.8)
  beeswarm(GSVA_results[i,]~factor(variable_of_interest), add=T,
           method="hex", corral=c("wrap"), pch=19, col=c(alpha("firebrick4", 0.2), alpha("gray50", 0.2)))
}
  dev.off()
}