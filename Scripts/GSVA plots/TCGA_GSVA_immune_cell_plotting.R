setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
library(beeswarm)
library(scales)

diseases <- list.files(recursive=F) #REMOVES STAD

for(x in diseases){
	mut_files<-	grep("cBio_MutsigCV_BB",list.files(paste(".",x,"Mutation",sep="/"),full.names=TRUE),value=TRUE)
	gsva_files<- grep("GSVA_RNASeq_log2.RData",list.files(paste(".",x,"Results",sep="/"),full.names=TRUE),value=TRUE)

	##Import and merge mutation files
	if(length(mut_files)==0){
		#print(paste(x,"does not have mutation data",sep=" "))
		next}

	if(length(gsva_files)==0){
		#print(paste(x,"does not have no GSVA Results",sep=" "))
		next}

	if(length(mut_files)==2){
		#print(paste(x,"has two mutation files",sep=" "))
		mutation_1 <- read.delim(mut_files[1],sep="\t",skip=1,row.names=1,header=TRUE)
		mutation_1 <- mutation_1[,-dim(mutation_1)[2]]
		mutation_2 <- read.delim(mut_files[2],sep="\t",skip=1,row.names=1,header=TRUE)
		mutation_2 <- mutation_2[,-dim(mutation_2)[2]]
		mutation<-cbind(mutation_1,mutation_2)
		mutation <- mutation[,-dim(mutation)[2]]
		}

	if(length(mut_files)==1){
		#print(paste(x,"has one mutation file",sep=" "))
		mutation <- read.delim(mut_files,sep="\t",skip=1,row.names=1,header=TRUE)
		mutation <- mutation[,-dim(mutation)[2]]
		}

	##Import and merge mutation and GSVA results files
	load(gsva_files)
	gsva_results<-gsva_results$es.obs
	colnames(gsva_results)<- apply(do.call(rbind,strsplit(colnames(gsva_results), split="-"))[,1:3],1,function(x){paste(x,collapse="-")})

	samples <-intersect(rownames(t(gsva_results)),rownames(mutation))
	mutation_set <- mutation[samples,]
	gsva_results_set <- t(gsva_results)[samples,]

	if(length(dim(mutation_set))!=2) {
		print(paste("There is a problem with the mutation data in dataset",x,sep=" "))
		next
	}

	if(length(dim(gsva_results_set))!=2) {
		print(paste("There is a problem with the GSVA data in dataset",x,sep=" "))
		next
	}

	if(dim(mutation_set)[1]!=dim(gsva_results_set)[1]) {
		print(paste("There is a general problem with Dataset",x,sep=" "))
		next
	}

	#print(paste("All is good with",x,sep=" "))
	system(paste("mkdir",paste(".",x,"Results","RNASeq_log2_gsva_immmune_cell",sep="/"),sep=" "))
	for(j in 1:dim(mutation_set)[2]){
		mut_interest<-mutation_set[,j]
		if(length(levels(factor(mut_interest)))==1){
		print(paste("No interesting mutation data for",names(mutation_set)[j],"in",x,sep=" "))
		next
		}
		if(length(levels(factor(mut_interest)))!=1){
		print(paste("Lots of interesting mutation data for",names(mutation_set)[j],"in",x,sep=" "))
		pdf(file = paste(paste(".",x,"Results","RNASeq_log2_gsva_immmune_cell",paste(colnames(mutation_set)[j],".pdf",sep=""),sep="/")))
		par(mfrow = c(2,2))
		mut_interest <- ifelse(mut_interest=="NaN","Wild-Type","Mutant")
		for(w in 1:length(colnames(gsva_results_set))){
			boxplot(gsva_results_set[,w]~factor(mut_interest),main=colnames(gsva_results_set)[w],
					xlab="Mutation Status",ylab="GSVA Relative Score",
					names=c("Mutant","Wild-Type"),cex.axis=0.9)
			beeswarm(gsva_results_set[,w]~factor(mut_interest),method="hex",  corral = c("omit"),
						col=c(alpha("red",0.5),
							alpha("black",0.5)),pch=19,
						add=TRUE)
		}
		dev.off()

		}

	}

}





