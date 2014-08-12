setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")
library(survival)
diseases <- list.files(recursive=F)[-c(10,15)] #REMOVES STAD
x <- "SKCM"
#for(x in diseases){
  y <- tolower(x)
  system(paste("mkdir", paste(".", x, "Plots", sep="/"), sep=" "))
  
  patient_clinical <- read.delim(paste(".", x, "CD", paste("nationwidechildrens.org_clinical_patient_", y,".txt", sep=""), sep="/"))
  patient_clinical <- patient_clinical[,c("bcr_patient_barcode", "ajcc_pathologic_tumor_stage", "vital_status","last_contact_days_to", "death_days_to")]
  rownames(patient_clinical) <- patient_clinical[,1]
  patient_clinical <- patient_clinical[-c(1:2),-1] #Removes patient barcode and 2 lines of heading
  patient_clinical[,"death_days_to"] <- as.numeric(gsub("[Not Applicable]",NA,patient_clinical[,"death_days_to"]))
  patient_clinical[,"last_contact_days_to"] <- as.numeric(gsub("[Not Available]",NA,patient_clinical[,"last_contact_days_to"]))

  days <- data.frame(alive = patient_clinical[,"last_contact_days_to"], dead = patient_clinical[,"death_days_to"])
  time <- cbind.data.frame(a=rownames(patient_clinical), time= rowSums(days, na.rm=T))
  time <- as.matrix(time[,-1])
  rownames(time) <- rownames(patient_clinical)

  patient_clinical <- cbind(time, patient_clinical)
  patient_clinical <- patient_clinical[,c("time", "ajcc_pathologic_tumor_stage", "vital_status")] # removes death days and last contact  
  writeLines(paste("Loaded Clinical", x, Sys.time()))
  
  load(paste(".", x, "Results", "LassoRegression_Mouse_BLZ945.RData", sep="/"), verbose=T)
  survival_exprs <- data.frame(merge(exprs.predict, patient_clinical, by="row.names"))
  rownames(survival_exprs) <- survival_exprs[,"Row.names"]
  survival_exprs <- survival_exprs[,-1]
  colnames(survival_exprs) <- c("Class","Score","Time", "Stage", "Status") 

  survival_exprs[,"Status"] <- ifelse(survival_exprs[,"Status"]=="Alive",0,
                                        ifelse(survival_exprs[,"Status"]=="Dead",1,NA))
  survival_exprs[,"Time"] <- as.numeric(gsub("[Not Applicable]",NA,survival_exprs[,"Time"]))
  survival_exprs[,"Time"] <- as.numeric(survival_exprs[,"Time"]/(365/12))
  survival_exprs[,"Class"] <- as.character(survival_exprs[,"Class"])
  survival_exprs[,"Score"] <- as.numeric(as.character(survival_exprs[,"Score"]))
  survival_exprs[,"Stage"] <- as.character(survival_exprs[,"Stage"])
  
  writeLines(paste("Survival Matrix Completed for", x, Sys.time()))
  
  types <- names(table(survival_exprs[,"Stage"])[])
  for(z in types){
      survival <- survfit(Surv(Time,Status)~Class, data=survival_exprs, subset=Stage==z)
      survmed <- as.matrix(summary(survival)$table[, "median"])
      rownames(survmed) <- c("BLZ945-like", "Vehicle-like")
      survmed <- survmed[1,]-survmed[2,]
      survmed <- round(as.numeric(survmed), digits=3)
      
      sdf <- survdiff(Surv(Time,Status)~Class, data=survival_exprs, subset=Stage==z)
      survp <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
      survp <- round(as.numeric(survp), digits=3)
      
      coxpvalue <- summary(coxph(data=survival_exprs,Surv(Time,Status)~Score, subset=Stage==z))$waldtest[3]
      coxpvalue <- round(as.numeric(coxpvalue), digits=3)
      
      pdf(file=paste(".", x, "Plots", paste("Stage_KM_",z, x, ".pdf", sep=""), sep="/"))
      plot(survival, main=paste(z, x, "Kaplan-Meier"),xlab="Months from Diagnosis to Death", ylab="Percent Overall Survival", col=c("firebrick3", "black"), lwd=2)
      legend("bottomleft", legend=c("Cox p-value:", coxpvalue,"Survival p-value:", survp,"Survival Difference BLZ945/Vehicle:",survmed), bty="n")
      dev.off()
    }
}