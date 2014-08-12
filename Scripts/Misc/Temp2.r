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

folder <- "ssGSEA_cyokine_CNA"
setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data/")

		#If so, load those files and format colnames
		CNA_mat <- read.delim(paste(".", x, "CNA", "all_thresholded.by_genes.txt", sep="/"))
		rownames(CNA_mat) <- CNA_mat[,1]
		CNA_mat <- CNA_mat[,-c(1:3)]
		colnames(CNA_mat) <- gsub("\\.","-",colnames(CNA_mat))
		colnames(gsva_results) <- gsub("\\.","-",colnames(gsva_results))

		colnames(CNA_mat) <- apply(do.call(rbind,strsplit(colnames(CNA_mat),split="-"))[,1:3],1,paste,sep="-",collapse="-")
		genes <- t(as.matrix(c("CXCL14","IL33","CCL27","IL26","IL25","TGFB3","CX3CL1","BMP7","TGFB1","CXCL14")))
		
		#Sample and genes selection
		genes <- intersect(rownames(CNA_mat),t(genes)[,1])
		samples <- intersect(colnames(CNA_mat),colnames(gsva_results))
		
		CNA_genes <- CNA_mat[genes,samples]
		gsva_results <- gsva_results[,samples]

		#make plots for each CNA event
		  writeLines(paste("Plotting for", x, Sys.time()))
		system(paste("mkdir",paste(".", x, "Plots", folder,sep="/")))
		  for(j in 1:length(genes)){ #loops over all mutations
		    mut_interest <- as.matrix(CNA_genes[j,])
		    
		    if(length(levels(factor(mut_interest)))==1){
		      print(paste("No interesting mutation data for",names(variable_of_interest)[j],"in",x,sep=" "))
		      next
		    }
		    if(length(levels(factor(mut_interest)))!=1){
		     
		        pdf(file=paste("./", x, "/Plots/", folder, "/", rownames(CNA_genes)[j],".pdf",sep=""))
		        par(mfrow=c(2,2))
		        for(i in 1:length(rownames(gsva_results))){ #loops over modules and cytokines
		              color_index <- round(seq.int(from=1,to=14,by=14/length(levels(factor(mut_interest)))))
							cols <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", 
						"#E6F598", "#ABDDA4", "#66C2A5", "aquamarine4","#3288BD", muted("dodgerblue2"), 
						 "#5E4FA2", muted("blueviolet"))[color_index]
		              boxplot(gsva_results[i,]~factor(mut_interest),
		                      main=rownames(gsva_results)[i],
		                      xlab="Variable of interest",ylab="Relative enrichment score",
		                      names=levels(factor(mut_interest)), cex.axis=0.9,las=0,
		                      col=cols)
		            }
		        }
		        dev.off()
		      }
		    
