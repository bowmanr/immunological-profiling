

setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
library(beeswarm)
library(scales)
library(caroline)
library(RColorBrewer)

library(ggplot2)
diseases <- list.files(recursive=F)#[c(1,10,11,15)]
x = "LGG"
folder <- "ssGSEA_immune_CNA"
for(x in diseases){
		#Determine if disease / folder has the appropriate files
		files <-	list.files(paste(".", x,"CNA",sep="/"))
		if(length(grep("GISTIC_genes.txt",files))!=1) {print(paste("No GISTIC genes file in",x)); next}
		  writeLines(paste("Loaded", x, sep=" "))
		files <-	list.files(paste("./", x,"/CNA",sep=""))
		if(length(grep("all_thresholded.by_genes.txt",files))!=1) {print(paste("No CNA file in",x)); next}
		  writeLines(paste("Loaded", x, sep=" "))
		files <-	list.files(paste(".", x,"Results",sep="/"))
		if(length(grep("RNASeq_counts_ssgsea_immmune_cell.RData",files))!=1) {print(paste("No ssGSEA results file in",x)); next}
		  writeLines(paste("Loaded", x, sep=" "))

		#If so, load those files and format colnames
		CNA_mat <- read.delim(paste(".", x, "CNA", "all_thresholded.by_genes.txt", sep="/"))
		rownames(CNA_mat) <- CNA_mat[,1]
		CNA_mat <- CNA_mat[,-c(1:3)]
		colnames(CNA_mat) <- gsub("\\.","-",colnames(CNA_mat))
		colnames(CNA_mat) <- apply(do.call(rbind,strsplit(colnames(CNA_mat),split="-"))[,1:3],1,paste,sep="-",collapse="-")
		load(paste(".", x, "Results", "RNASeq_counts_ssgsea_immmune_cell.RData",sep="/"),verbose=TRUE)
		colnames(gsva_results) <- apply(do.call(rbind,strsplit(colnames(gsva_results),split="-"))[,1:3],1,paste,sep="-",collapse="-")		
		genes <- read.delim(paste(".", x,"CNA","GISTIC_genes.txt",sep="/"),sep=" ",header=FALSE)

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
		                      col=cols)}
		        }
		        dev.off()
		      }
		    }














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
