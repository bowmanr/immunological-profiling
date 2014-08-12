library(matrixStats)
library(beeswarm)
library(caroline)
library(scales)
library(gdata)
setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data/")
load("/Users/bowmanr/Projects/Imm_profile/Datasets/Gene_sets/cytokines_of_interest.RData",verbose=T)

		#Determine if disease / folder has the appropriate files
		files <-	list.files(paste("./", x,"/CNA",sep=""))
		if(length(grep("all_thresholded.by_genes.txt",files))!=1) {print(paste("No CNA file in",x)); next}
		  writeLines(paste("Loaded", x, sep=" "))
		files <-	list.files(paste("./", x,"/RNASeq",sep=""))
		if(length(grep("RNASeqCounts.RData",files))!=1) {print(paste("No RNASeq file in",x)); next}
		  writeLines(paste("Loaded", x, sep=" "))

		#If so, load those files and format colnames
		CNA_mat <- read.delim(paste(".", x, "CNA", "all_thresholded.by_genes.txt", sep="/"))
		rownames(CNA_mat) <- CNA_mat[,1]
		CNA_mat <- CNA_mat[,-c(1:3)]
		colnames(CNA_mat) <- apply(do.call(rbind,strsplit(colnames(CNA_mat),split="\\."))[,1:3],1,paste,sep=".",collapse=".")

		setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
		setwd("/Users/bowmanr/Projects/BLZ945/Analysis/12.1.12")
		load("./exprs.list.RData",verbose=TRUE)
		RNASeqFinal <- exprs.list[[i]]
		setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")

		#Sample and cytokine selection
		genes <- Reduce(intersect,list(rownames(CNA_mat),rownames(RNASeqFinal),cytokines_of_interest))
		samples <- intersect(colnames(CNA_mat),colnames(RNASeqFinal))
		
		#Subset CNA matrix on cytokines and samples of interest
		CNA_cytokine <- CNA_mat[genes,samples]
		RNA_cytokine <- RNASeqFinal[genes,samples]

		## Merge the CNA and RNA data
		merged <- list()
		for(i in 1:length(rownames(CNA_cytokine))) {
			merged[[i]]	<-t(rbind(CNA_cytokine[i,],RNA_cytokine[i,]))
			colnames(merged[[i]]) <- c("CNA","RNA")
		}
		names(merged) <- rownames(CNA_cytokine)

		## Plot the results

		pdf(file=paste("./", x, "/Plots/", "/cytokines_microarray_box.pdf",sep=""))
		par(mfrow=c(2,2))

		for(i in 1:length(merged)){
			color_index <- round(seq.int(from=1,to=14,by=14/length(levels(factor(merged[[i]][,"CNA"])))))
			cols <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "aquamarine4","#3288BD", muted("dodgerblue2"), "#5E4FA2", muted("blueviolet"))[color_index]

		  	boxplot(merged[[i]][,"RNA"]~factor(merged[[i]][,"CNA"]), main=names(merged)[i], xlab="Copy Number Alteration (GISTIC)",ylab="Relative Expression (log2)", col=cols)
		 	#violins(merged[[i]][,"RNA"], by=levels(factor(merged[[i]][,"CNA"])), main=names(merged)[i], xlab="Copy Number Alteration (GISTIC)",ylab="RNASeq Counts", deciles=FALSE, connect="median", connectcol="black",drawRect=T, col=cols)
			  #axis(1, labels = FALSE)
			  #text(1:length(color_index), par("usr")[3]-0.03, srt = 45, adj = 1, labels = levels(factor(variable_of_interest)), xpd=T, cex=0.9) 
			  #mtext(1, text = "Stage", line = 3.5)
		}
		dev.off()
