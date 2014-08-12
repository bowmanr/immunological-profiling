setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data")

geneList <- read.delim("./ANOVA-baseline vs. other 28 conditions.txt")
symbols <- geneList[,1:2]
meanRem <- grep("Mean", colnames(geneList), value=T, invert=T)
geneList <- geneList[meanRem]
rownames(geneList) <- geneList[,"Probeset.ID"]
geneList <- geneList[,-1]

names <- gsub("Fold.Change.baseline.vs..","Cyto=",colnames(geneList))
names <- grep("Cyto=", names, value=T)

is.even <- function(x) x%%2==0
P_FClist <- list()

counter <- 1
for(x in seq(from=2, to=length(colnames(geneList)), by=2)){
  P_FClist[[counter]] <- geneList[,x:(x+1)]
  counter <- counter+1
}
names(P_FClist) <- names

pvalcut <- 0.1
fccut <- 2


Results <- lapply(P_FClist,function(x){
              cbind(x,p.adjust(x[,1],method="fdr",n=length(x[,1])))
            })

Results_test <- lapply(Results,function(x){
  t(apply(x, 1, function(y){
      c(as.numeric(y[3]) < pvalcut,
            y[2] > fccut,
            y[2] < fccut*(-1))
    }))
})



resIndex <- lapply(Results_test, function(x){
  t(apply(x, 1, function(y){
    cbind(y[1]&y[2], y[1]&y[3])
  }))
})

names(resIndex) <- names

resIndex <- lapply(resIndex,function(x){
  colnames(x) <- c("Up","Down")
  return(x)
})
      
res_symbols <- lapply(resIndex, function(x){
    merge(x,symbols,by.x="row.names",by.y="Probeset.ID")
})


res_set <- lapply(res_symbols,function(x){
  list("up"=c(x[x[,4]!=""&x[,2]==TRUE,4]),"down"=c(x[x[,4]!=""&x[,3]==TRUE,4]))
})

do.call(rbind,lapply(res_set,function(x){list(length(x[[1]]),length(x[[2]]))}))


save(res_set,file="CytoList_final_gene_sets.RData")