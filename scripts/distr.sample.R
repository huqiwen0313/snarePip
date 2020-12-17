# This script is used to merge samples from different libraries based on ATAC/RNA folder structure
# Usage: Rscript merge.samples.R sample_dir TRUE/FALSE RNA/ATAC TRUE/FALSE snare_2 path_to_sample_table
# Qiwen Hu - 2020

library(SnapATAC)
library(pagoda2)
library(snarePip)

##
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
# delete tmp file or not
tmpFile <- args[2]
# RNA or ATAC
type <- args[3]
# if link sample ID to libID
linkLIbID <- args[4]
# assay type
assay <- args[5]
# path to sample table
link.path <- args[6]
samtool.path <- args[7]

tmp.path <- file.path(path, "tmp")
if(type=="ATAC"){
  directory <- file.path(tmp.path, "snap")
} else{
  directory <- file.path(tmp.path, "pagoda_RData")
}

# directory of link table
if(linkLIbID){
  link.table <- read.table(link.path, sep="\t", header=TRUE)
}

files <- unique(gsub("[.].*|-SP.*|_N.*", "", dir(directory)))
for(i in 1:length(files)){
  tissue <- unique(snarePip:::getTissue(link.table, files[i]))
  filePath <- file.path(path, assay, tissue, "samples")
  if(!dir.exists(file.path(path, assay, tissue))){
    dir.create(file.path(path, assay))
    dir.create(file.path(path, assay, tissue))
    dir.create(file.path(path, assay, tissue, "samples"))
  }
  
  # sample snap files
  if(type=="ATAC"){
    snapFiles <- system(paste("ls", file.path(directory, paste(files[i], "*.snap", sep="")), sep=" "), 
                        intern=TRUE)
    # extract fragment file paths
    frag.path <- file.path(tmp.path, "fragments")
    fragmentFiles <- system(paste("ls", file.path(frag.path, paste(files[i], "*.sort.bed", sep="")), sep=" "), 
                            intern=TRUE)
  
    # create folders
    system(paste("mkdir", file.path(filePath, files[i]), sep=" "),
           intern=TRUE)
    system(paste("mkdir", file.path(filePath, files[i], "Sample_output"), sep=" "),
           intern=TRUE)
    system(paste("mkdir", file.path(filePath, files[i], "Sample_output", "snap"), sep=" "),
           intern=TRUE)
    system(paste("mkdir", file.path(filePath, files[i], "Sample_output", "objects"), sep=" "),
           intern=TRUE)
    system(paste("mkdir", file.path(filePath, files[i], "Sample_output", "fragments"), sep=" "),
           intern=TRUE)
    system(paste("mkdir", file.path(filePath, files[i], "Sample_output", "macs2"), sep=" "),
           intern=TRUE)
    system(paste("mkdir", file.path(filePath, files[i], "Sample_output", "bam"), sep=" "),
           intern=TRUE)
    
    # merge fragment files 
    system2(command="cat", args=c(paste(fragmentFiles, collapse = ' '),
                                  ">", paste(file.path(filePath, files[i], "Sample_output", "fragments", files[i]), 
                                             "fragments.bed", sep=".")))
    # extract bam file paths
    bam.path <- file.path(tmp.path, "bam")
    bamFiles <- system(paste("ls", file.path(bam.path, paste(files[i], "*.bam", sep="")), sep=" "), 
                            intern=TRUE)
    bamFiles <- bamFiles[grep(paste(files[i], ".*rmsk.bam", sep=""), bamFiles, perl=TRUE)]
    bamFiles <- paste(bamFiles, collapse=" ")
   
    # merge bamFiles
    outputbamFile <- paste(file.path(filePath, files[i], "Sample_output", "bam", files[i]), "bam", sep=".")
    system(paste(samtool.path, "merge", outputbamFile, bamFiles, "-@ 15", sep=" "),
           intern=TRUE)
    
    # create initial merged p2 object
    p2 <- Pagoda2$new()
    p2[["fragmentsPaths"]] <- paste(file.path(filePath, files[i], "Sample_output", "fragments", files[i]), 
                                    "fragments.bed", sep=".")
    saveRDS(p2, paste(file.path(filePath, files[i], "Sample_output", "objects", files[i]), "p2.obj.rds", sep="."))
  } else{
    pagodaFiles <- system(paste("ls", file.path(directory, paste(files[i], "*filtered.RData", sep="")), sep=" "), 
                          intern=TRUE)
    system(paste("mkdir", file.path(filePath, files[i]), sep=" "),
           intern=TRUE)
    system(paste("mkdir", file.path(filePath, files[i], "Sample_output"), sep=" "),
           intern=TRUE)
    system(paste("mkdir", file.path(filePath, files[i], "Sample_output", "obj"), sep=" "),
           intern=TRUE)
    cMatrix <- lapply(pagodaFiles, function(r){
      t <- load(r)
      get(t)})
    mergedMatrix <- snarePip:::merge.sparse(cMatrix)
    saveRDS(mergedMatrix, paste(file.path(filePath, files[i], "Sample_output", "obj", files[i]), "sample_matrix.rds", sep="."))
  }
  
  # Split experiment data
  expDir <- file.path(file.path(filePath, files[i]), "Experiment_output")
  system(paste("mkdir", expDir, sep=" "), intern=TRUE)
  if(type=="ATAC"){
    folders <- c("bam", "fastqFiles", "fragments", "macs", "objects",
                 "pmats", "QCs", "reports", "snap")
  } else{
    folders <- c("alignment", "dropest_out", "FastQC", "pagoda_RData", 
                 "QCs", "seurat_obj", "tagged", "reports")
  }
  
  for(j in 1:length(folders)){
    system(paste("mkdir", file.path(expDir, folders[j]), sep=" "), intern=TRUE)
    if(linkLIbID){
      #folder.files <- dir(file.path(tmp.path, folders[j]))
      folder.files <- system(paste("ls", file.path(tmp.path, folders[j], paste(files[i], "*", sep="")), sep=" "),
                             intern=TRUE)
      folder.files <- gsub(".*/", "", folder.files)
      
      lapply(1:length(folder.files), function(r){
        # get full name experiment ID
        sampleid <- gsub("_S\\d+.*", "_S", folder.files[r])
        sampleid <- gsub("[.].*", "", folder.files[r])
        libID <- snarePip:::getLibID(link.table, sampleid)
        if(length(libID) > 0){
          surfix <- gsub(".*_S\\d+.", "", folder.files[r])
          surfix <- gsub(".*\\d+.", "", folder.files[r])
          sample.files <- file.path(tmp.path, paste(folders[j], folder.files[r], sep="/"))
          target.files <- file.path(expDir, 
                                    paste(folders[j], paste(libID, surfix, sep="."), sep="/"))
          system(paste("cp", sample.files, target.files, sep=" "))
          #if(folders[j] == "QCs"){
          #  target.files <- file.path(expDir, 
          #                            paste(folders[j], folder.files[r], sep="/"))
          #  system(paste("cp", sample.files, target.files, sep=" "))
          #}
        } else{
          sample.files <- file.path(tmp.path, paste(folders[j], folder.files[r], sep="/"))
          target.files <- file.path(expDir, 
                                    paste(folders[j], folder.files[r], sep="/"))
          system(paste("cp", sample.files, target.files, sep=" "))
        }
      })} else{
      system(paste("cp", paste(file.path(tmp.path, folders[j]), "/", files[i], "*", sep=""),
                   file.path(expDir, folders[j]), sep=" "), intern=TRUE)
    }
  }
}

# move log files out
system(paste("mv", file.path(tmp.path, "logs/*"), file.path(path, "logs"), sep=" "), intern=TRUE)
if(tmpFile){
  system(paste("rm -rf", file.path(path, "tmp"), sep=" "), intern=TRUE)
}
system(paste("rm", file.path(path, "*filtered.bam"), sep=" "), intern=TRUE)
