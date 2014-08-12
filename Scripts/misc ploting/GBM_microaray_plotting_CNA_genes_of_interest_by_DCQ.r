setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
library(GSEABase)
library(GSVAdata) 
library(GSVA)
library(parallel)
library(marray)
library(ggplot2)
library(scales)

#Xue Cytokine List
load("/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/misc data/CytoList_final_gene_sets.RData", verbose=T)
finalset <- unlist(res_set, recursive=F)
#Galon-Engler
load("/Users/bowmanr/Projects/Imm_profile/Datasets/Gene_sets/Galon_Engler_list.RData", verbose=T)

mod_set <- galon_engler
genes <- unique(c("CX3CL1","CXCL14","CXCL12","CSF1","BMP8A","BMP8B","CXCL16","CCL2","IL15","TNFSF13B","TGFB1","BMP1","CCL28","BMP4","TGFB3","TNFSF11","IL15","TNFSF10","TGFBI"))
diseases <- list.files(recursive=F)#[c(1,10,11,15)]
folder <- "ssGSEA_Cytokines_interest_lung_CNA"

for(x in c("LUAD","LUSC")){
		files <-	list.files(paste("./", x,"/CNA",sep=""))
		if(length(grep("all_thresholded.by_genes.txt",files))!=1) {print(paste("No CNA file in",x)); next}
		files <-	list.files(paste("./", x,"/RNASeq",sep=""))
		if(length(grep("RNASeqCounts.RData",files))!=1) {print(paste("No RNASeq file in",x)); next}
		
		load(paste(".", x, "RNASeq", "RNASeqCounts.RData", sep="/"), verbose=T)
		writeLines(paste("Loaded", x, sep=" "))

		gsva_results <- gsva(expr=as.matrix(RNASeqFinal), gset.idx.list=mod_set, method="ssgsea", rnaseq=TRUE,
		                       min.sz=0, max.sz=10000, parallel.sz = 12, verbose=TRUE)
		colnames(gsva_results) <- gsub("-",".",colnames(gsva_results))

		CNA_mat <- read.delim(paste(".", x, "CNA", "all_thresholded.by_genes.txt", sep="/"))
		rownames(CNA_mat) <- CNA_mat[,1]
		CNA_mat <- CNA_mat[,-c(1:3)]
		colnames(CNA_mat) <- apply(do.call(rbind,strsplit(colnames(CNA_mat),split="\\."))[,1:3],1,paste,sep=".",collapse=".")
			
		samples <- intersect(colnames(gsva_results),colnames(CNA_mat))		
		CNA_genes <- CNA_mat[genes,samples]
		gsva_results <- gsva_results[,samples]

		#make plots for each CNA event
		writeLines(paste("Plotting for", x, Sys.time()))
		system(paste("mkdir",paste(".", x, "Plots", folder,sep="/")))
		  for(j in 1:length(genes)){ #loops over all mutations
		    mut_interest <- as.matrix(CNA_genes[j,])
		    
		    if(length(levels(factor(mut_interest)))==1){
		      print(paste("No interesting mutation data for",colnames(CNA_genes)[j],"in",x,sep=" "))
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
		   }
		    
