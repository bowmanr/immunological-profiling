######################
### Libraries used
######################
source("/Users/bowmanr/Projects/Imm_profile/Scripts/DCQ/Immune_cell_abundance_functions.r")

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


load("/Volumes/MacintoshHD3/data/TCGA_Assembler/Diseases_non_TCGA/Skin/GSE14905/GSE14905pheno.RData",verbose=TRUE)
load("/Volumes/MacintoshHD3/data/TCGA_Assembler/Diseases_non_TCGA/Skin/GSE14905/GSE14905exprs.RData",verbose=TRUE)

  	gsva_results <- gsva(expr=as.matrix(exprs_1), gset.idx.list=galon_engler, method="ssgsea", rnaseq=FALSE,
                       min.sz=0, max.sz=10000, parallel.sz = 24, verbose=TRUE)

	gsva_results<-regression_input(exprs=exprs_1,
											immune_exprs=Allantaz,
											test_gene_set_random=FALSE,
											limit_tumor_genes=FALSE,
											test_gene_set=Allantaz_genes,
											alphaParam=0.05,
											quant_norm=FALSE,
											lambda_param=0.1,
											method="elastic")
#mut_interest <- factor(pheno_1[,"Tissue"],levels=c("normal","pilocytic astrocytoma","ependymoma","glioblastoma","medulloblastoma"))
mut_interest <- factor(pheno_1[,"Disease"],levels=c("Normal","Uninvolved","Lesion"))

color_index <- round(seq.int(from=1,to=14,by=14/length(levels(factor(mut_interest)))))
cols <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "aquamarine4","#3288BD", muted("dodgerblue2"), "#5E4FA2", muted("blueviolet"))[color_index]


keyword <- "elastic_allantaz"
pdf(file=paste("/Users/bowmanr/Desktop/Final Plots/Temp/",keyword,".pdf",sep=""))
for(i in 1:length(rownames(gsva_results))){ #loops over modules and cytokines
      boxplot(gsva_results[i,]~factor(mut_interest),
              main=rownames(gsva_results)[i],
              xlab="Variable of interest",ylab="Relative enrichment score",
              names=levels(factor(mut_interest)), cex.axis=0.9,las=0,
              col=cols)
  		axis(1, labels = FALSE)
		text(1:length(color_index), par("usr")[3]-0.03, srt = 45, adj = 1, labels = levels(factor(mut_interest)), xpd=T, cex=0.9) 
		mtext(1, text = "Stage", line = 3.5)}
      dev.off()