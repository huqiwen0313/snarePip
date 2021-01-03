# this script is used to generate sample-level QC object and statistics
# Qiwen Hu - 2020

library(pagoda2)
library(dropestr)
library(Seurat)
library("DropletUtils")
suppressMessages(library(org.Hs.eg.db))
library(snarePip)

## functions
GetGenesetFraction_overall <- function (count.matrix, genes){
  umi.counts <- sort(Matrix::colSums(count.matrix), decreasing = T)
  presented.mit.genes <- intersect(genes, rownames(count.matrix))
  genes.frac <- sum(Matrix::colSums(count.matrix[presented.mit.genes, 
                                                 names(umi.counts)]))/sum(umi.counts)
  return(genes.frac)
}
####

# parameters
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
sampleID <- args[2]
scrublet_path <- args[3]

# QC cutoffs
min.cell.size <- 200
min.ngenes <- 200
max.ngenes <- 7500
min.ncells <- 500

objdir <- file.path(path, "obj")
count <- readRDS(file.path(objdir, paste(sampleID, "sample_matrix.rds", sep=".")))
# filtering
sizePcell <- colSums(count)
# filtering by cell size
fcells <- names(sizePcell[sizePcell>min.cell.size])
count <- count[, colnames(count) %in% fcells]


# doublet score
if(dim(count)[2] > min.cell.size){
  scrublescore <- snarePip:::GetScrubletScores(count, min.molecules.per.gene=1, method="scrublet",
                                    pythonPath=scrublet_path)
  doublet.detect <- snarePip:::GetScrubletScores(count, min.molecules.per.gene=10, method="doubletDetection",
                                      pythonPath=scrublet_path)
  if(!file.exists(file.path(path, "doublets"))){
    dir.create(file.path(path, "doublets"))
  }
  doublet.path <- file.path(path, "doublets")
  saveRDS(list(scrublet=scrublescore, doubletDetection=doublet.detect), 
          file.path(doublet.path, paste(sampleID, "doublet.rds", sep=".")))
}


# pagoda object
p2 <- snarePip:::p2proc(count, n.cores=20, min.cells.per.gene=5, min.cell.size=200, n.odgenes=2e3, 
                         get.largevis=FALSE, make.geneknn=FALSE, perplexity=50, log.scale = TRUE, 
                         batch = NULL, nPcs = 30, k = 30, get.tsne = TRUE)
if(dim(count)[2] > min.cell.size){
  p2[["doublet"]] <- list(scrublet=scrublescore, doubletDetection=doublet.detect)
}
saveRDS(p2, file.path(objdir, paste(sampleID, "p2.rds", sep=".")))

# create p2 App
if(dim(count)[2] > min.cell.size){
  ids <- unlist(lapply(mget(colnames(p2$counts), org.Hs.egALIAS2EG, ifnotfound=NA), function(x) x[1]))
  rids <- names(ids); names(rids) <- ids
  go.env <- list2env(eapply(org.Hs.egGO2ALLEGS,function(x) as.character(na.omit(rids[x]))))
  p2$makeGeneKnnGraph(n.cores=5) 
  p2web <- snarePip:::GetPagodaWebApp(p2, p2$clusters$PCA$multilevel, additional.metadata=list(), 
                verbose=T, go.sets=NULL, go.env=go.env, test.pathways=TRUE)
  # save object
  p2web$serializeToStaticFast(paste(objdir,"/",sampleID, ".p2App.bin",sep=""), verbose=T)
}

# generate QC statistics
mit_genes <- c("MT-ND1",
               "MT-ND2",
               "MT-ND3",
               "MT-ND4",
               "MT-ND4L",
               "MT-ND5",
               "MT-ND6",
               "MT-CYB",
               "MT-CO1",
               "MT-CO2",
               "MT-CO3",
               "MT-ATP6",
               "MT-ATP8")

aligned_UMIs_per_cell <- Matrix::colSums(count)
median_UMIs_per_cell <- median(aligned_UMIs_per_cell)
number_of_cells_after_filtering <- ncol(count)

if(number_of_cells_after_filtering>1){
  percent_read_mapping_to_mitochrondrial<- round(GetGenesetFraction_overall(count, mit_genes), 4)
} else{
  percent_read_mapping_to_mitochrondrial<-"Too few cells to calculate"
}

if(number_of_cells_after_filtering>1){
  median_fraction_read_in_mitochondrial <- round(median(GetGenesetFraction(count, mit_genes)), 4)
}else{
  median_fraction_read_in_mitochondrial <-"Too few cells to calculate"
}

if(number_of_cells_after_filtering < min.ncells){
  Sample_quality <- "Not_Pass"
} else{Sample_quality <- "Pass"}

qc.stat <- c(median_UMIs_per_cell, median_fraction_read_in_mitochondrial,
             percent_read_mapping_to_mitochrondrial, number_of_cells_after_filtering,
             Sample_quality, paste(objdir,"/",sampleID, ".p2App.bin",sep=""))

measurement <- c("median_UMIs_per_cell", "median_fraction_read_in_mitochondrial",
                 "percent_read_mapping_to_mitochrondrial", "number_of_cells_after_filtering",
                 "Sample_quality", "p2_webApp_address")


qc.stat <- data.frame(measurement=measurement, stat=qc.stat)
qc.dir <- file.path(path, "QCs")
if(!file.exists(qc.dir)){
  dir.create(qc.dir)
}
write.table(qc.stat, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t", 
            file.path(qc.dir, paste(sampleID, "qc.txt", sep=".")))

# output tenX format
tenX.dir <- file.path(path, "matrix_tenX")
write10xCounts(tenX.dir, count, overwrite=T)
