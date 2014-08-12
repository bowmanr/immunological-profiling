setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
library(matrixStats)
library(beeswarm)
library(caroline)
library(scales)
library(gdata)
###########################################
#Plotting all Gene Expression Based on Stage
###########################################
folder <- "Stage_ssGSEA_Macrophage" #Name of folder for results
data <- "ssGSEA_Cyto_Module_RNASeqCounts.RData" #Name of data input
diseases <- list.files(recursive=F)[-c(10,15)] #Removes STAD and Merged
for(x in diseases){
  y <- tolower(x)
  system(paste("mkdir", paste(".", x, "Plots", folder, sep="/"), sep=" "))
  
  gsva_files<- grep(data, list.files(paste(".",x,"Results", "ssGSEA_Cyto", sep="/"),full.names=TRUE),value=TRUE)
  load(gsva_files, verbose=T)
  
  stage <- read.delim(paste(".", x, "CD", paste("nationwidechildrens.org_clinical_patient_", y,".txt", sep=""), sep="/"))
  if(c("ajcc_pathologic_tumor_stage") %in% colnames(stage)){
  stage <- stage[,c("bcr_patient_barcode", "ajcc_pathologic_tumor_stage")]
  stage <- stage[-c(1:2),] #Removes 2 lines of heading
  rownames(stage) <- stage[,1]
  colnames(stage) <- c("Patient", "Stage")
  stageFinal <- as.matrix(stage[, "Stage"])
  stageFinal <- gsub("[ABC]", "", stageFinal[,1])
  
  if(x == "LUSC"){
    stageFinal1 <- as.matrix(stageFinal[stageFinal=="Stage I"|stageFinal=="Stage II"|stageFinal=="Stage III"])
    rownames(stageFinal1) <- rownames(stage)[stageFinal=="Stage I"|stageFinal=="Stage II"|stageFinal=="Stage III"]
    
    }else{
      stageFinal1 <- as.matrix(stageFinal[stageFinal=="Stage I"|stageFinal=="Stage II"|stageFinal=="Stage III"|stageFinal=="Stage IV"])
      rownames(stageFinal1) <- rownames(stage)[stageFinal=="Stage I"|stageFinal=="Stage II"|stageFinal=="Stage III"|stageFinal=="Stage IV"]      
  } #Not enough samples in Stage IV for voilins plot
############## Sets up Variables for Plotting
  index <- intersect(colnames(gsva_results), rownames(stageFinal1))
  
  variable_of_interest <- as.matrix(stageFinal1[index,])
  colnames(variable_of_interest) <- "Stage"
  
  GSVA_results <- gsva_results[,index]
  
  writeLines(paste("Plotting for", x, Sys.time()))
  
  color_index <- round(seq.int(from=1,to=14,by=14/length(levels(factor(variable_of_interest)))))
  cols <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "aquamarine4","#3288BD", muted("dodgerblue2"), "#5E4FA2", muted("blueviolet"))[color_index]
  
  #pvalue <- p.adjust(apply(GSVA_results,1,function(x)t.test(x~factor(variable_of_interest))$p.value), n=length(rownames(GSVA_results)))
  
  #if(any(pvalue<=(0.01))){
  
    for(j in 1:length(colnames(variable_of_interest))){
      pdf(file=paste("./", x, "/Plots/", folder, "/", colnames(variable_of_interest)[j],".pdf",sep=""))
      par(mfrow=c(2,2))
      
      for(i in 1:length(rownames(GSVA_results))){ #this loops over each module
        #if(pvalue[i]<=(0.01)){ 
            violins(GSVA_results[i,], by=factor(variable_of_interest[,j]), main=rownames(GSVA_results)[i], xlab="",ylab="Relative enrichment score", deciles=FALSE, connect="median", connectcol="black", names=rep(c(""),length(color_index)),drawRect=T, col=cols)
            axis(1, labels = FALSE)
            text(1:length(color_index), par("usr")[3]-0.02, srt = 45, adj = 1, labels = levels(factor(variable_of_interest)), xpd=T, cex=0.9) 
            mtext(1, text = "Stage", line = 3.5)
            #}
      }
      dev.off()
    }
  } else{next}
}

