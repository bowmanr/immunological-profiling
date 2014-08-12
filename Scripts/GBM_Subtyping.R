## Necessary Packages:
library(kernlab)
library(glmnet)
library(matrixStats)
#Import data sets
setwd("/Users/bowmanr/Projects/publicDataSets")
TCGA.clinical <- read.delim("./rawData/clinical/TCGA.clinical.patient.gbm.txt") # limited clinical info
TCGA.clacnc.terms <- as.matrix(read.table("./rawData/misc/CLACNC.terms.txt")) # List of 840 genes used to define the 4 subtypes
TCGA.core <- read.delim("./rawData/misc/TCGA.core.txt") # This matrix contain info whether samples were used to define the core subtpes in Verhaak et al 2008

TCGA.core[,"Sample"] <- gsub("\\.", "-", TCGA.core[,"Sample"])

load("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data/GBM/RNASeq/zScore.RData", verbose=T)
TCGAnorm <- exprs_zScore

# Identify samples with clinical information and expression information
# Then subset those matrices based on these samples.
rownames(TCGA.core) <- TCGA.core[,"Sample"]

coreSamples <- TCGA.core[TCGA.core[,"Core"]=="CORE",]

## Subset the initial TCGA matrix to only contain "CORE" samples
#Check for NA's
## Mean cetner this matrix such that all gens are on the same scale.  I do not standarize with standard deviation.
trainVector <- coreSamples[, "Subtype"]

table(rownames(coreSamples)%in%colnames(TCGAnorm))

trainMatrix <- TCGA.core[]
# Set the rest of the samples in a matrix where we will learn their subtypes

## This block is used to subset the clinical info such that only core samples are listed.

## A factor is then made of their subtpes, such that we now have 3 pieces of information:
## 1) We have a matrix to train on: (TCGA.train)
## 2) We have a factor with the training information (Training)
## 3) We have a matrix to learn the subtypes on: (TCGA.learn)
## This is then fed into a support vector machine (SVM) using the "ksvm" function from the "kernlab" package

## SVM  
svm.TCGA<-ksvm(t(TCGA.train), #Training matrix
               Training, # "Truth" factor
               cross=10, # 10 fold cross validation to learn which SVM is best
               kernel="vanilladot",  # vanilladot kernel keeps things simple
               family="multinomial", # we are predicting 4 classes
               prob.model=TRUE,  # We want a probability model included so we can give each predicted sample a score, not just a class
               scale=FALSE)  # We already scaled before by mean-centering the data.

svm.TCGA  # Output

## Run the SVM on the prediction dataset, and tidy it into a nice matrix
TCGA.subtype.call      <- as.matrix(predict(svm.TCGA, t(TCGA.learn)))
TCGA.subtype.prob      <- as.matrix(predict(svm.TCGA, t(TCGA.learn ), type="probabilities"))
TCGA.subtype.all       <- data.frame(TCGA.subtype.call,TCGA.subtype.prob)
rownames(TCGA.subtype.all) <- rownames(t(TCGA.learn))
TCGA.subtype.all[,"Sample"] <- rownames(TCGA.subtype.all)

## write to file
write.table(TCGA.subtype.all,"TCGA.subtype.non.core.txt",sep="\t")

## OTher matrix merging.  I can't remember why I did all of this, sorry.
TCGA.core.samples.2 <- TCGA.core.samples.1[,c(1,2)]
TCGA.clinical.subset <- TCGA.clinical[,c(2,9,31)]
TCGA.subtype.subset<- TCGA.subtype.all[,c("Sample","TCGA.subtype.call")] ; colnames(TCGA.subtype.subset)<-c("Sample","Subtype")
TCGA.subtype.complete <- rbind(TCGA.core.samples.2,TCGA.subtype.subset)

TCGA.clinical.int <- merge(TCGA.subtype.complete,TCGA.clinical.subset,by.x="row.names",by.y="row.names")
TCGA.clinical.int[,"status"] <- ifelse(TCGA.clinical.int[,"vital_status"]=="DECEASED",1,0)
TCGA.clinical.final <- TCGA.clinical.int[,c(2,5,7,3)]; colnames(TCGA.clinical.final)<-c("Sample","survival","status","Subtype")

setwd("/Users/bowmanr/Projects/publicDataSets/rawData/clinical")
write.table(TCGA.clinical.final,"TCGA.clinical.data.1.txt",sep="\t")