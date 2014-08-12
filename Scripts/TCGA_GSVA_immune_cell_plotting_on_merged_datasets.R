setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
library(beeswarm)
library(scales)
library(caroline)
library(RColorBrewer)
library(ggplot2)

	load("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data/Merged/Results_ssgsea_counts_immune_All_RNASeqCounts.RData",verbose=TRUE)
	load("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data/Merged/All_samples_pheno.RData",verbose=TRUE)

	rownames(All_samples) <- All_samples[,1]
	samples <-intersect(rownames(t(gsva_results)),rownames(All_samples))
	All_samples <- All_samples[samples,]
	gsva_results_set <- t(gsva_results)[samples,]
	mut_interest <- droplevels(factor(All_samples[,2]))

color_index <- c(1:14)
cols <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "aquamarine4","#3288BD", muted("dodgerblue2"), "#5E4FA2", muted("blueviolet"))[color_index]

		pdf(file = paste("./Merged/Results/RNASeq_counts_immune_All.pdf",sep=""))
		par(mfrow = c(2,1))
		for(w in 1:length(colnames(gsva_results_set))){
			violins(gsva_results_set[,w],by=factor(mut_interest), main=colnames(gsva_results_set)[w],
			xlab="", ylab="ssGSEA Score", deciles=FALSE, connect="median", connectcol="black", col=cols, names=rep(c(""),length(color_index)),drawRect=T)
			axis(1, labels = FALSE)
      text(1:length(color_index), par("usr")[3]-0.02, srt = 45, adj = 1, labels = levels(mut_interest), xpd=T, cex=0.9) 
			mtext(1, text = "Mutation Status", line = 3.5)
			
				}
		dev.off()
