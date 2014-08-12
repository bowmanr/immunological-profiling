setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")


##################### Start Loop over Diseases ############################
diseases <- list.files(recursive=F) #REMOVES STAD and Merged
for(x in diseases){
    results_files<-  grep("Results.RData",list.files(paste(".",x,"Results","GSEA",sep="/"),full.names=TRUE),value=TRUE)

	if(length(results_files)==0){
		print(paste("No results for",x,"move on",sep=" "))
		next
	}

	load(results_files,verbose=TRUE)
	cor_mat<-data.frame(do.call(cbind,lapply(gsea_results_final,function(x){
	ifelse(x=="NULL",c(rep("NULL",105)),x[,"Correlation"])})))

	direction_mat <-data.frame(do.call(cbind,lapply(gsea_results_final,function(x){
	ifelse(x=="NULL",c(rep("NULL",105)),x[,"Direction"])})))

	pval_mat <-data.frame(do.call(cbind,lapply(gsea_results_final,function(x){
	ifelse(x=="NULL",c(rep("NULL",105)),x[,"FDR"])})))

}