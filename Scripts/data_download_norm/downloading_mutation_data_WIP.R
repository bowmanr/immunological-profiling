

################ Loading TCGA ###############
setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler")

source("./Scripts/Module_A.r")
source("./Scripts/Module_B.r")
################ Clinical Data ################
master_folder <- "/Volumes/MacintoshHD3/data/TCGA_Assembler/Data"
setwd("/Volumes/MacintoshHD3/data/TCGA_Assembler/Data/")
load("/Volumes/MacintoshHD3/data/TCGA_Assembler/Scripts/DirectoryTraverseResult_Jan-30-2014.rda")

diseases <- list.files(recursive=F)
for(cancerType in diseases){
  # load directory traverse result and create folder to store output files
  writeLines("Identifying MAF files.");
  disease_IDs =  grep(pattern = toupper(cancerType), x = upper_file_url, ignore.case = FALSE,value=TRUE);
  SpecificID  =  grep(pattern = toupper(paste("mutations", sep = "")), x = disease_IDs, ignore.case = FALSE,value=TRUE);
  maf_files   =  grep(pattern =".MAF$", x = SpecificID, ignore.case = FALSE,value=TRUE);
  indices<-match(maf_files,upper_file_url)
  # search for the Sample and Data Relationship Format (SDRF) file of the specified platform and cancer type
  if (length(indices) == 0)
  {
    writeLines("No MAF files detected"); 
    next();
  }
  if (length(indices) > 1)
  {
   writeLines("Downloading information from identified MAF file.");

    URL = GetNewestURL(AllURL = file_url[indices]);
  }
  downloadResult = urlReadTable(url = URL);
  if (downloadResult$errorFlag != 0)
  {
    writeLines("Error in downloading SDRF file.");
    next();
  }
institute_index <- grep(".edu",downloadResult$data[1:400,])
Hugo_Symbol_index <- grep("Hugo_Symbol",downloadResult$data[1:400,])
if((institute_index[1]-institute_index[2])==(institute_index[1]-institute_index[2])){
  start <- institute_index[1]-2
  stop <- institute_index[2]-2
  ncols <- (institute_index[2]-institute_index[1])
  col_name_set <- downloadResult$data[c(Hugo_Symbol_index:(institute_index[1]-3)),]
if(cancerType=="GBM"){
  ncols=78
}
  mutation_info = matrix(toupper(downloadResult$data[start:dim(downloadResult$data)[1],1]),ncol=ncols,byrow=T)
  
  colnames(mutation_info) <- col_name_set[1:dim(mutation_info)[2]]
}else{
  writeLines("This file is weird, get out while you can");
  next();
}
system(paste("mkdir ",paste(master_folder,cancerType,"Mutation",sep="/")))
save(mutation_info,file=paste(master_folder,cancerType,"Mutation","Mutation_assembler.RData",sep="/"))
   writeLines(paste("Completed", cancerType,sep=" "));
}



