######################
### Libraries used
######################
library(Hmisc)
library(survival)
library(matrixStats)
#####################
##  Data Import
#####################
options(stringsAsFactors=FALSE)
setwd("/Users/bowmanr/Projects/BLZ945/Analysis/12.1.12")
load("./exprs.list.RData")
load("./clinical.list.RData")
setwd("/Users/bowmanr/Projects/publicDataSets/")
load("./Analysis/Classification/1.26.13.TAM.lasso.sig.RData",verbose=T)
load("./Analysis/Classification/6.19.13.TAM.SVM.sig.RData",verbose=T)


GCIMP<-read.delim("/Users/bowmanr/Projects/publicDataSets/rawData/misc/TCGA.GCIMP.txt", header=T)
GCIMP.Class<-merge(list.exprs.lasso[[1]],GCIMP,by.x="row.names",by.y="Sample")
length(which(GCIMP.Class$GCIMP=="GCIMP.1"))
(which(GCIMP.Class$GCIMP=="GCIMP.1"&!GCIMP.Class$Class=="BLZ"))
Holder <- as.list(NULL)
for(i in 1:5){
	exprs.1 <- merge(clinical.list[[i]], list.exprs.lasso[[i]], by.x="Sample", by.y="row.names")
	#rownames(exprs.1) <- exprs.1[,1]
	Holder[[i]] <- exprs.1
}
names(Holder) <- names(clinical.list)
Genes <- c("ADM")
Subtypes <- c("Proneural","Mesenchymal","Neural","Classical")
subtype <- c("Proneural")

Final <- list(NULL,NULL)
names(Final) <- c("TCGA","Combination")
for( i in 1:5){
 	samples <- as.character((Holder[[i]])[,"Sample"])#[which((Holder[[i]])[,"Subtype"]==subtype)]
	exprs <- (exprs.list[[i]])[Genes,samples]
	exprs <- (exprs -rowMeans(exprs))
	exprs.t <- t(exprs)
	exprs.1.t <- merge(Holder[[i]],exprs.t, by.x="Sample", by.y="row.names")
	exprs.1.t[,"Score"] <- as.numeric(as.character(exprs.1.t[,"Score"]))
	exprs.1.t[,"BLZ.Score"]<- 1-as.numeric(as.character(exprs.1.t[,"Score"]))
	if(i==1){
		Final[[1]] <-exprs.1.t
			}
	else{
		Final[[2]] <- rbind(Final[[2]],exprs.1.t)
	}
}

Final[[1]][,"survival.m"] <- Final[[1]][,"survival"]/(365/12)
Final[[2]][,"survival.m"] <- Final[[2]][,"survival"]/(365/12)

Fit<-survfit(Surv(survival.m,status)~Class,  data=Final[[1]])
Fit.p <- 1 - pchisq(survdiff(Surv(survival.m,status)~Class, data=Final[[1]])$chisq,1)

#pdf(file = paste("/Volumes/MacintoshHD3/data/TCGA_Assembler/Plots for Poster/Finished Plots", "Whole_GBM_KM.pdf",sep="/"))
par(#mfrow=c(2,2),
  oma=c(2,2,2,2))
plot(Fit,lwd=3,xlab="Months",ylab="Percent Survival",
     col=c("red","black"), lty=c(1),cex.axis=1.5,cex.lab=1.5,
     main=c("Lasso prediction \n in TCGA-proneural"),cex.main=1.5,font.main=1)
legend("topright", c("Vehicle-Like","BLZ945-Like"), cex=1.4, 
       bty="n", lty=1,col=c("black","red","white"),lwd=3)

#dev.off()










Fit.proneural<-survfit(Surv(survival.m,status)~Class,  subset=Subtype=="Proneural", data=Final[[1]])
Fit.proneural.p <- 1 - pchisq(survdiff(Surv(survival.m,status)~Class,  subset=Subtype=="Proneural", data=Final[[1]])$chisq,1)

Fit.Mesenchymal<-survfit(Surv(survival.m,status)~Class,  subset=Subtype=="Mesenchymal", data=Final[[1]])
Fit.Mesenchymal.p <- 1 - pchisq(survdiff(Surv(survival.m,status)~Class,  subset=Subtype=="Mesenchymal", data=Final[[1]])$chisq,1)

Fit.Classical<-survfit(Surv(survival.m,status)~Class,  subset=Subtype=="Classical", data=Final[[1]])
Fit.Classical.p <- 1 - pchisq(survdiff(Surv(survival.m,status)~Class,  subset=Subtype=="Classical", data=Final[[1]])$chisq,1)

Fit.Neural<-survfit(Surv(survival.m,status)~Class,  subset=Subtype=="Neural", data=Final[[1]])
Fit.Neural.p <- 1 - pchisq(survdiff(Surv(survival.m,status)~Class,  subset=Subtype=="Neural", data=Final[[1]])$chisq,1)

par(#mfrow=c(2,2),
  oma=c(2,2,2,2))
plot(Fit.proneural,lwd=3,xlab="Months",ylab="Percent Survival",
     col=c("red","black"), lty=c(3,1),cex.axis=1.5,cex.lab=1.5,
     main=c("Lasso prediction \n in TCGA-proneural"),cex.main=1.5,font.main=1)
legend("topright", c("Vehicle-Like","BLZ945-Like",paste("p <",signif(Fit.proneural.p,3))), cex=1.4, 
       bty="n", lty=c(1,3),col=c("black","red","white"),lwd=3)



plot(Fit.Mesenchymal,lwd=3,xlab="Months",ylab="Percent Survival",
     col=c("red","black"), lty=c(3,1),cex.axis=1.5,cex.lab=1.5,
          main=c("Lasso prediction \n in TCGA-mesenchymal"),cex.main=1.5,font.main=1)
legend("topright", c("Vehicle-Like","BLZ945-Like",paste("p <",signif(Fit.Mesenchymal.p,3))), cex=1.4,  
       bty="n", lty=c(1,3),col=c("black","red","white"),lwd=3)


plot(Fit.Classical,lwd=3,xlab="Months",ylab="Percent Survival",
     col=c("red","black"), lty=c(3,1),cex.axis=1.5,cex.lab=1.5,
          main=c("Lasso prediction \n in TCGA-classical"),cex.main=1.5,font.main=1)
legend("topright", c("Vehicle-Like","BLZ945-Like",paste("p <",signif(Fit.Classical.p,3))), cex=1.4, 
       bty="n", lty=c(1,3),col=c("black","red","white"),lwd=3)

plot(Fit.Neural,lwd=3,xlab="Months",ylab="Percent Survival",
     col=c("red","black"), lty=c(3,1),cex.axis=1.5,cex.lab=1.5,
          main=c("Lasso prediction \n in TCGA-neural"),cex.main=1.5,font.main=1)
legend("topright", c("Vehicle-Like","BLZ945-Like",paste("p <",signif(Fit.Neural.p,3))), cex=1.4, 
       bty="n", lty=c(1,3),col=c("black","red","white"),lwd=3)







