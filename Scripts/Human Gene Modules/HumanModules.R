library(gdata)

setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts")

roughGenes <- read.xls("./Human Gene Modules/TableS2B.xlsx", sheet="Table S2B", verbose=F, method="csv", blank.lines.skip=T)

justMods <- roughGenes[-c(1,2),-1]

modSep <- apply(justMods, MARGIN=2, FUN=list)
mods <- lapply(modSep, function(x){
               unique(x[[1]])})
names(mods) <- c(paste("Module", 1:49, sep="_"))
noBlank <- lapply(mods, function(x){
                  x[x!=""]})

save(noBlank, file = paste("./Human Gene Modules/ModuleLists.RData"))