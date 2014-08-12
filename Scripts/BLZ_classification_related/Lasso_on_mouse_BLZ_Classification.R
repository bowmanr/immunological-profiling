######################
### Libraries used
######################
library(Hmisc)
library(survival)
library(coin)
library(matrixStats)
library(beeswarm)
library(ggplot2)

#####################
##  Human Data Import
#####################
setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
load("/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/Human_Gene_Lists.RData", verbose=T)

################################
##  Mouse Data Import and setup
################################
mouse_exprs <- read.delim("/Users/bowmanr/Projects/BLZ945/Final/CR03.exprs.data.txt",sep="\t")

#Gene Lists
human_to_mouse <- as.matrix(read.delim("/Users/bowmanr/Projects/generalScripts/human_to_mouse.txt", header=T))
mouse_genes <- as.matrix(rownames(mouse_exprs))
human_genes <- as.matrix(Gene.List)

rownames(human_to_mouse)<- as.matrix(human_to_mouse[,"Mouse"])

mouseHomol <- as.matrix(human_to_mouse[human_to_mouse[,"Mouse"]%in%mouse_genes,])
allHomol <- as.matrix(mouseHomol[mouseHomol[,"Human"]%in%human_genes,])

Mouse.1 <- mouse_exprs[allHomol[,"Mouse"],]
Mouse.1 <-Mouse.1[!is.na(Mouse.1[,1]),]
Mouse.1 <-(Mouse.1 - rowMeans(Mouse.1)) 
Mouse.1.t <- t(Mouse.1)

################################
## Lasso Training
################################
library(glmnet)

target <- as.factor(c(rep("Vehicle",8),rep("BLZ945",8)))

cv.glmnet<-cv.glmnet(Mouse.1.t,target, family="binomial",  standardize=FALSE, nfolds=4)
glmnet.1<-glmnet(Mouse.1.t,target, family="binomial",  standardize=FALSE, lambda=c(cv.glmnet$lambda.min))

A<-as.matrix(glmnet.1$beta)
names(A) <- colnames(Mouse.1.t)
B<- as.matrix(A)[which(A[,1]!=0.0)]

################################
## Disease Loading and Formatting
################################
diseases <- list.files(recursive=F)[-14] #REMOVES STAD
for(x in diseases){
  y <- tolower(x)
  
  load(paste(".", x, "RNASeq", "zScore.RData", sep="/"), verbose=T) #DNE for STAD
  exprs_data <- t(exprs_zScore[allHomol[,"Human"],])
  
  drug_clinical <- read.delim(paste(".", x, "CD", paste("nationwidechildrens.org_clinical_drug_", y,".txt", sep=""), sep="/"))
  patient_clinical <- read.delim(paste(".", x, "CD", paste("nationwidechildrens.org_clinical_patient_", y,".txt", sep=""), sep="/"))
  
  writeLines(paste("Loaded", x))
  
  ################################
  ## Lasso Prediction
  ################################
  exprs.class <- predict(glmnet.1, newx=exprs_data, type="class")
  exprs.score <- predict(glmnet.1, newx=exprs_data, type="response")
  rownames(exprs.class) <- rownames(exprs_data)
  rownames(exprs.score) <- rownames(exprs_data)
  exprs.predict <- cbind(exprs.class,exprs.score)
  colnames(exprs.predict) <- c("Class","Score")
  
  save(exprs.predict,file=paste(".", x, "Results", "LassoRegression_Mouse_BLZ945.RData", sep="/"))
  print(table(exprs.class))
}

#Gene.List <- rownames(exprs_zScore)
#save(Gene.List, file="./Scripts/Human_Gene_Lists.RData")
#Used to create the common human genes