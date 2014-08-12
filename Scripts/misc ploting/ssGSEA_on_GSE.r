################GSVA on Normal, non-TCGA datasets ##########################
setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Diseases_non_TCGA/GSVA_Results")
library(beeswarm)
library(scales)
library(caroline)
library(RColorBrewer)
library(ggplot2)

pheno_wd <- "/Volumes/MacintoshHD3/data/TCGA_Assembler/Diseases_non_TCGA"

GSE_files <- grep("GSE", list.files(), value=T)
for(x in GSE_files){
	    GSE_number <-   do.call(rbind,strsplit(x,split="_"))[,2]
	    GSE_disease <-  do.call(rbind,strsplit(x,split="_")) [,1]
	    GSE_keyword <-  do.call(rbind,strsplit(x,split="_")) [,3]
	    load(x, verbose=T)
	    gsva_results_set <- gsva_results
	    pheno_files <- grep("pheno", list.files(paste(pheno_wd,GSE_disease,GSE_number,sep="/"),full.names=TRUE), value=T)
	    load(pheno_files, verbose=T)

	    for(j in 2:length(colnames(pheno_1))){
	    		mut_interest <- factor(pheno_1[,j])

				color_index <- round(seq.int(from=1,to=14,by=14/length(levels(factor(mut_interest)))))
				cols <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", 
						"#E6F598", "#ABDDA4", "#66C2A5", "aquamarine4","#3288BD", muted("dodgerblue2"), 
						 "#5E4FA2", muted("blueviolet"))[color_index]

	    		pdf(file = paste(pheno_wd,GSE_disease,GSE_number,paste("ssGSEA_results_",GSE_keyword,"_",colnames(pheno_1)[j],".pdf",sep=""),sep="/"))
				par(mfrow = c(2,2))
					for(w in 1:length(rownames(gsva_results_set))){
						boxplot(gsva_results_set[w,]~factor(mut_interest), main=rownames(gsva_results_set)[w], 
							xlab="",names=F,
							ylab="ssGSEA Enrichment Score", col=cols)
					  axis(1, labels = FALSE)
					  text(1:length(color_index), par("usr")[3]-0.005, srt = 45, adj = 1, labels = levels(factor(mut_interest)), xpd=T, cex=0.9) 
					 # mtext(1, text = "Stage", line = 3.5)
					}
					dev.off()
					}
	}