setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Diseases_non_TCGA")
source("/Users/bowmanr/Projects/generalScripts/affy_utils_1.R")

##############################################
####BRCA GSE
##############################################
disfold <- "./Breast/"
x <- "GSE7904"
pheno <- read.delim(file=paste(disfold,x,"/GSE7904_series_matrix.txt", sep=""), skip = 29, nrows = 66-31, header=T)
pheno[,2]
pheno_1 <- t(pheno[c(1,7,9),])
rownames(pheno_1) <- pheno_1[,1]
pheno_1 <- pheno_1[-1,]
colnames(pheno_1) <- c("Sample_CEL","Sample_ID","Group")

data <- read.delim(file=paste(disfold,x,"/GSE7904_series_matrix.txt", sep=""), skip = 66, header=T, row.names=1)
exprs <- log2(data)
exprs_1<-summarize.exprs.data(exprs,rownames="hgu133plus2.db")

save(exprs_1, file=paste(disfold, x, "exprs.RData", sep=""))
save(pheno_1, file=paste(disfold, x, "pheno.RData", sep=""))
##############################################
x <- "GSE15852"
pheno <- read.delim(file=paste(disfold,x,"/GSE15852_series_matrix.txt", sep=""), skip = 37, nrows = 79-38, header=T)
pheno[,2]
pheno_1 <- t(pheno[c(1,7,10,11),])
rownames(pheno_1) <- pheno_1[,1]
pheno_1 <- pheno_1[-1,]
pheno_1[,3] <- do.call(rbind,strsplit(pheno_1[,3],split="histopathological exam: "))[,2]
pheno_1[,4] <- do.call(rbind,strsplit(pheno_1[,4],split="grade: "))[,2]

colnames(pheno_1) <- c("Sample_CEL","Tissue_source","Exam","Grade")

data <- read.delim(file=paste(disfold,x,"/GSE15852_series_matrix.txt", sep=""), skip = 80, header=T, row.names=1)
exprs <- log2(data)
exprs_1<-summarize.exprs.data(exprs,rownames="hgu133a.db")

save(exprs_1, file=paste(disfold, x, "exprs.RData", sep=""))
save(pheno_1, file=paste(disfold, x, "pheno.RData", sep=""))
##############################################
####LUNG GSE
##############################################
disfold <- "./Lung/"
x <- "GSE2514"
pheno <- read.delim(file=paste(disfold,x,"/GSE2514-GPL8300_series_matrix.txt", sep=""), skip = 40, nrows = 69-40, header=T)
pheno[,2]
pheno_1 <- t(pheno[c(1,11),])
rownames(pheno_1) <- pheno_1[,1]
pheno_1 <- pheno_1[-1,]
colnames(pheno_1) <- c("Sample_CEL","Tissue")

data <- read.delim(file=paste(disfold,x,"/GSE2514-GPL8300_series_matrix.txt", sep=""), skip = 69, header=T, row.names=1)
exprs <- log2(data)
exprs_1<-summarize.exprs.data(exprs,rownames="hgu95av2.db") #GPL8300 = hgu95av2.db WILL NOT RUN

save(exprs_1, file=paste(disfold, x, "exprs.RData", sep=""))
save(pheno_1, file=paste(disfold, x, "pheno.RData", sep=""))
##############################################
x <- "GSE10245"
pheno <- read.delim(file=paste(disfold,x,"/GSE10245_series_matrix.txt", sep=""), skip = 35, nrows = 67-35, header=T)
pheno[,2]
pheno_1 <- t(pheno[c(1,9),])
rownames(pheno_1) <- pheno_1[,1]
pheno_1 <- pheno_1[-1,]
colnames(pheno_1) <- c("Sample_CEL","Tissue")
pheno_1[,2] <- do.call(rbind,strsplit(pheno_1[,2],split="disease state: "))[,2]

data <- read.delim(file=paste(disfold,x,"/GSE10245_series_matrix.txt", sep=""), skip = 67, header=T, row.names=1)
exprs <- data
exprs_1<-summarize.exprs.data(exprs,rownames="hgu133plus2.db")

save(exprs_1, file=paste(disfold, x, "exprs.RData", sep=""))
save(pheno_1, file=paste(disfold, x, "pheno.RData", sep=""))
##############################################
####BRAIN GSE
##############################################
disfold <- "./Brain/"
x <- "GSE15824"
pheno <- read.delim(file=paste(disfold,x,"/GSE15824_series_matrix.txt", sep=""), skip = 35, nrows = 75-35, header=T)
pheno[,2]
pheno_1 <- t(pheno[c(1,9,10,13),])
rownames(pheno_1) <- pheno_1[,1]
pheno_1 <- pheno_1[-1,]
colnames(pheno_1) <- c("Sample_CEL","Tissue","Group","Grade")
pheno_1[grep("cell line:",pheno_1[,2]),4] <- c(": Cell line")
pheno_1[grep("astrocytes",pheno_1[,3]),4] <- c(": Astrocytes")
pheno_1[pheno_1[,3]=="",3] <- c(": Normal brain")
pheno_1[grep("Normal brain",pheno_1[,3]),4] <- c(": Normal brain")
pheno_1[grep("cell line:",pheno_1[,2]),3] <- c(": Cell line")
pheno_1[grep("cell line:",pheno_1[,2]),4] <- pheno_1[grep("cell line:",pheno_1[,2]),3]
pheno_1[,2] <- do.call(rbind,strsplit(pheno_1[,2],split=": "))[,2]
pheno_1[,3] <- do.call(rbind,strsplit(pheno_1[,3],split=": "))[,2]
pheno_1[,4] <- do.call(rbind,strsplit(pheno_1[,4],split=": "))[,2]


data <- read.delim(file=paste(disfold,x,"/GSE15824_series_matrix.txt", sep=""), skip = 75, header=T, row.names=1)
exprs <- log2(data)
exprs_1<-summarize.exprs.data(exprs,rownames="hgu133plus2.db")

save(exprs_1, file=paste(disfold, x, "exprs.RData", sep=""))
save(pheno_1, file=paste(disfold, x, "pheno.RData", sep=""))
##############################################
x <- "GSE32374"
pheno <- read.delim(file=paste(disfold,x,"/GSE32374_series_matrix.txt", sep=""), skip = 35, nrows = 72-35, header=T)
pheno[,2]
pheno_1 <- t(pheno[c(37,10),])
rownames(pheno_1) <- pheno_1[,1]
pheno_1 <- pheno_1[-1,]
colnames(pheno_1) <- c("Sample_CEL","Tissue")
pheno_1[,2] <- do.call(rbind,strsplit(pheno_1[,2],split=": "))[,2]

data <- read.delim(file=paste(disfold,x,"/GSE32374_series_matrix.txt", sep=""), skip = 72, header=T, row.names=1)
exprs <- data
exprs_1<-summarize.exprs.data(exprs,rownames="hgu133plus2.db")

save(exprs_1, file=paste(disfold, x, "exprs.RData", sep=""))
save(pheno_1, file=paste(disfold, x, "pheno.RData", sep=""))
##############################################
x <- "GSE50161"
pheno <- read.delim(file=paste(disfold,x,"/GSE50161_series_matrix.txt", sep=""), skip = 29, nrows = 65-29, header=T)
pheno[,2]
pheno_1 <- t(pheno[c(1,9,10),])
rownames(pheno_1) <- pheno_1[,1]
pheno_1 <- pheno_1[-1,]
colnames(pheno_1) <- c("Sample_CEL","Disease","Tissue")
pheno_1[,2] <- do.call(rbind,strsplit(pheno_1[,2],split=": "))[,2]
pheno_1[,3] <- do.call(rbind,strsplit(pheno_1[,3],split=": "))[,2]
pheno_1[grep("normal",pheno_1[,2]),3] <- pheno_1[grep("normal",pheno_1[,2]),2]

data <- read.delim(file=paste(disfold,x,"/GSE50161_series_matrix.txt", sep=""), skip = 66, header=T, row.names=1)
exprs <- data
exprs_1<-summarize.exprs.data(exprs,rownames="hgu133plus2.db")

save(exprs_1, file=paste(disfold, x, "exprs.RData", sep=""))
save(pheno_1, file=paste(disfold, x, "pheno.RData", sep=""))

##############################################
####Melanoma GSE
##############################################
disfold <- "./Skin/"
x <- "GSE46517"
pheno <- read.delim(file=paste(disfold,x,"/GSE46517_series_matrix.txt", sep=""), skip = 72, nrows = 93-72, header=T)
pheno[,2]
pheno_1 <- t(pheno[c(21,7),])
rownames(pheno_1) <- pheno_1[,1]
pheno_1 <- pheno_1[-1,]
colnames(pheno_1) <- c("Sample_CEL","Disease")
pheno_1[,2] <- do.call(rbind,strsplit(pheno_1[,2],split="\\."))[,1]
pheno_1[,2] <- gsub("Met","Metastatic",pheno_1[,2] )
pheno_1[,2] <- gsub("Pri","Primary",pheno_1[,2] )
pheno_1[,2] <- gsub("Nev","Benign Nevus",pheno_1[,2] )
pheno_1[,2] <- gsub("Cont","Control Skin",pheno_1[,2] )

data <- read.delim(file=paste(disfold,x,"/GSE46517_series_matrix.txt", sep=""), skip = 93, header=T, row.names=1)
exprs <- log2(data)
exprs_1<-summarize.exprs.data(exprs,rownames="hgu133a.db")

save(exprs_1, file=paste(disfold, x, "exprs.RData", sep=""))
save(pheno_1, file=paste(disfold, x, "pheno.RData", sep=""))


##############################################
####Melanoma GSE
##############################################
disfold <- "./Skin/"
x <- "GSE14905"
pheno <- read.delim(file=paste(disfold,x,"/GSE14905_series_matrix.txt", sep=""), skip = 29, nrows = 65-29, header=T)
pheno[,2]
pheno_1 <- t(pheno[c(1),])
pheno_1 <- cbind(pheno_1,rownames(pheno_1))
rownames(pheno_1) <- pheno_1[,1]
pheno_1 <- pheno_1[-1,]
colnames(pheno_1) <- c("Sample_CEL","Disease")
pheno_1[,2] <- do.call(rbind,strsplit(pheno_1[,2],split="\\."))[,1]

data <- read.delim(file=paste(disfold,x,"/GSE14905_series_matrix.txt", sep=""), skip = 65, header=T, row.names=1)
exprs_1<-summarize.exprs.data(data,rownames="hgu133plus2.db")

save(exprs_1, file=paste(disfold, x, "exprs.RData", sep=""))
save(pheno_1, file=paste(disfold, x, "pheno.RData", sep=""))
