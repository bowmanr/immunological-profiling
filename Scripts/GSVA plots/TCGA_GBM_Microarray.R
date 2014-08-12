######################
### Libraries used
######################
library(Hmisc)
library(survival)
library(matrixStats)
library(GSEABase)
library(GSVAdata) 
library(GSVA)
library(parallel)
library(marray)
library(ggplot2)
library(caroline)
library(scales)
#Xue Modules
load("/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/Human Gene Modules/ModuleLists.RData", verbose=T)
#Galon-Engler
load("/Users/bowmanr/Projects/Imm_profile/Datasets/Gene_sets/Galon_Engler_list.RData", verbose=T)
#Xue Cytokine List
load("/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/misc data/CytoList_final_gene_sets.RData", verbose=T)
finalset <- unlist(res_set, recursive=F)
#Merged Cytokine and Mod List
mergedlist <- galon_engler

#####################
##  Data Import
#####################
options(stringsAsFactors=FALSE)
setwd("/Users/bowmanr/Projects/BLZ945/Analysis/12.1.12")
load("./exprs.list.RData",verbose=TRUE)
load("./clinical.list.RData",verbose=TRUE)


modlist <- mergedlist #See above for options

  gsva_results <- gsva(expr=as.matrix(exprs.list[[1]]), gset.idx.list=modlist, method="ssgsea", rnaseq=FALSE,
                       min.sz=0, max.sz=10000, parallel.sz = 24, verbose=TRUE)


	rownames(clinical.list[[1]]) <- clinical.list[[1]][,1]
	samples <-intersect(rownames(t(gsva_results)),clinical.list[[1]][,1])
	gsva_results_set <- gsva_results[,samples]
	clinical.list[[1]] <- clinical.list[[1]][samples,]

	mut_interest <-	factor(clinical.list[[1]][,"Subtype"])

		              color_index <- round(seq.int(from=1,to=14,by=14/length(levels(factor(mut_interest)))))
							cols <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", 
						"#E6F598", "#ABDDA4", "#66C2A5", "aquamarine4","#3288BD", muted("dodgerblue2"), 
						 "#5E4FA2", muted("blueviolet"))[color_index]	;
						 	setwd("/Users/bowmanr/Desktop")
		pdf(file = paste("test.pdf",sep=""))
		par(mfrow = c(2,2))
		for(w in 1:length(rownames(gsva_results_set))){
		              boxplot(gsva_results_set[w,]~factor(mut_interest),
		                      main=rownames(gsva_results)[w],
		                      xlab="Variable of interest",ylab="Relative enrichment score",
		                      names=levels(factor(mut_interest)), cex.axis=0.9,las=0,
		                      col=cols)}
			
		dev.off()
}
