source("/Users/bowmanr/Projects/Imm_profile/Scripts/DCQ/Immune_cell_abundance_functions.r")
setwd("/Users/bowmanr/Projects/BLZ945/Analysis/12.1.12")
load("./exprs.list.RData",verbose=TRUE)
load("./clinical.list.RData",verbose=TRUE)
setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")

x="GBM"
i=1
folder <- "elastic_allantaz_mutation"
system(paste("mkdir", paste(".", x, "Plots", folder, sep="/"), sep=" "))
  
  mut_files<-  grep("cBio_MutsigCV_BB",list.files(paste(".",x,"Mutation",sep="/"),full.names=TRUE),value=TRUE)
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




setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
		exprs <- exprs.list[[i]]
		pheno <- clinical.list[[i]][,c("Sample","Subtype")]
		rownames(pheno) <- pheno[,1]
		samples <- intersect(rownames(pheno),colnames(exprs))
		exprs <- exprs[,samples]
		pheno_1 <- as.matrix(pheno[samples,-1])
		rownames(pheno_1) <- pheno[,"Sample"]

		regression_results<-regression_input(exprs=exprs,
											immune_exprs=Allantaz,
											test_gene_set_random=FALSE,
											limit_tumor_genes=FALSE,
											test_gene_set=Allantaz_genes,
											alphaParam=0.05,
											quant_norm=FALSE,
											lambda_param=0.1,
											method="elastic")

		gsva_results <- regression_results
		colnames(gsva_results) <- gsub("\\.","-",colnames(gsva_results))
	index <- intersect(colnames(gsva_results), rownames(mutation))
  
  mutation_set <- mutation[index,]
  GSVA_results <- gsva_results[,index]
  
  writeLines(paste("Plotting for", x, Sys.time()))
  
  for(j in 1:length(colnames(mutation_set))){ #loops over all mutations
    mut_interest<-mutation_set[,j] 
    mut_interest <- ifelse(mut_interest=="NaN","Wild-Type","Mutant")
    
    if(length(levels(factor(mut_interest)))==1){
      #print(paste("No interesting mutation data for",names(variable_of_interest)[j],"in",x,sep=" "))
      next
    }
    if(length(levels(factor(mut_interest)))!=1){
      #print(paste("Lots of interesting mutation data for",names(variable_of_interest)[j],"in",x,sep=" "))
      
      if(any(as.matrix(table(mut_interest))[,1]<3)){
        next
      }
      
           
        pdf(file=paste("./", x, "/Plots/", folder, "/", colnames(mutation_set)[j],".pdf",sep=""))
        par(mfrow=c(2,2))
        for(i in 1:length(rownames(GSVA_results))){ #loops over modules and cytokines
                               color_index <- round(seq.int(from=1,to=14,by=14/length(levels(factor(mut_interest)))))
              cols <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", 
            "#E6F598", "#ABDDA4", "#66C2A5", "aquamarine4","#3288BD", muted("dodgerblue2"), 
             "#5E4FA2", muted("blueviolet"))[color_index]
                  boxplot(GSVA_results[i,]~factor(mut_interest),
                          main=rownames(GSVA_results)[i],
                          xlab="Variable of interest",ylab="Relative enrichment score",
                          names=levels(factor(mut_interest)), cex.axis=0.9,las=0,
                          col=cols) 
        }
        dev.off()
      }}