# re-call peaks at sample-level
library(pagoda2)
library(snarePip)
library(GenomicRanges)
library(IRanges)

## arguments
args = commandArgs(trailingOnly=TRUE)
# path for ATAC processed data
ATACpath <- args[1]
# path for RNA processed data -  processed directory: ../RNA_processed_data
RNApath <- args[2]
# sampleID - data to be analyzed
sampleID <- args[3]
# path to macs
path.to.macs <- args[4]
# gsize e.g. hs, mm
gsize <- args[5]
assay <- args[6]

## read ref files
#link.dir <- "/d0/data/ucsd/refs"
link.dir <- args[7]
link.table <- read.table(args[7], sep="\t", header=TRUE)

tissue <- snarePip:::getTissue(link.table, sampleName=sampleID)
# read ATAC initial p2 object
atacP2file <- paste(file.path(ATACpath, assay, tissue, "samples", sampleID), "Sample_output", "objects",
                   paste(sampleID, "p2.obj.rds", sep="."), sep="/")
atacP2obj <- readRDS(atacP2file)
rnaCountfile <- paste(file.path(RNApath, assay, tissue, "samples", sampleID), "Sample_output", "obj",
                        paste(sampleID, "sample_matrix.rds", sep="."), sep="/")
rnaCount <- readRDS(rnaCountfile)
# generate p2 object
if(ncol(rnaCount) > 500){
  rnaP2obj <- snarePip:::p2proc(rnaCount)
  atacP2obj[["RNAp2obj"]] <- rnaP2obj
  saveRDS(atacP2obj, file=atacP2file)
} 

  
# call population based peaks
# temporary based on multilevel clusters
if(ncol(rnaCount) > 500){
  RNAclusters <- rnaP2obj$clusters$PCA$multilevel
} else{
  RNAclusters <- setNames(rep("1", length(rownames(atacP2obj[["pmat"]]))), 
                          rownames(atacP2obj[["pmat"]]))
}

out.dir <- paste(file.path(ATACpath, assay, tissue, "samples", sampleID), "Sample_output", "macs2", sep="/")
peaks.gr <- callMACsCluster(atacP2obj, clusters=RNAclusters, path.to.macs, gsize,
                              out.dir=out.dir,
                              macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
                              ncores=1, keep.simple=TRUE,
                              min.cells=10)





