# script for generate filtered count matrix
# Qiwen Hu & Xin Wang - 2020

require('httr')
library('broom')
require('igraph')
require('png')
require('rlang')
require('crayon')
require('digest')
require('assertthat')
require('glue')
require('purrr')
require('backports')
require('Seurat')
require('dplyr')
require('tidytext')
library("stringr")
library("magrittr")
library("dropestr")
library(pagoda2)
library(snarePip)

#related functions
merge.dtn6 <- function(rds.file, sampleID, dt.barcodes, n6.barcodes){
  rds <- readRDS(rds.file)
  tidy.cm <- tidy(rds$cm)  

  n6.dt.map <- dt.barcodes
  names(n6.dt.map) <- n6.barcodes
  
  dtn6.barcodes <- sapply(tidy.cm$column, function(x) substr(x, 17, 24))
  r3r2.barcodes <- sapply(tidy.cm$column, function(x) substr(x, 1, 16))

  is.dt <- sapply(dtn6.barcodes, function(x) x %in% dt.barcodes) 
  merged.barcodes <- sapply(dtn6.barcodes, function(x) {
                          if (x %in% names(n6.dt.map)) { 
                            return(n6.dt.map[[x]])
                          } else { return(x) }
                        })
  #add libID
  LibID <- snarePip:::getLibID(link.table, sampleID)
  new.cell.barcodes <- paste0(LibID,"_", r3r2.barcodes, merged.barcodes)
  tidy.cm$column <- new.cell.barcodes
  combined.matrix <- tidy.cm %>% group_by(row, column) %>% 
    dplyr::summarise(mergedValue = sum(value, na.rm = TRUE)) %>% cast_sparse(row, column, mergedValue)
  return(combined.matrix)
}


merge.dtn6.aligned.reads.umis.per.cell <- function(rds.file, sampleID, dt.barcodes, n6.barcodes){
  rds <- readRDS(rds.file)
  info <- data.frame(aligned_reads_per_cell=rds$aligned_reads_per_cell, 
                     aligned_umis_per_cell=rds$aligned_umis_per_cell)
  info$row <- rownames(info)
  n6.dt.map <- dt.barcodes
  names(n6.dt.map) <- n6.barcodes
  dtn6.barcodes <- sapply(rownames(info), function(x) substr(x, 17, 24))
  r3r2.barcodes <- sapply(rownames(info), function(x) substr(x, 1, 16))
  is.dt <- sapply(dtn6.barcodes, function(x) x %in% dt.barcodes) 
  merged.barcodes <- sapply(dtn6.barcodes, 
                            function(x) {
                              if (x %in% names(n6.dt.map)) { 
                                return(n6.dt.map[[x]])
                              } else { return(x) }
                            })
  
  #add libID
  LibID <- snarePip:::getLibID(link.table, sampleID)
  new.cell.barcodes <- paste0(LibID,"_",r3r2.barcodes, merged.barcodes)
  info$row <- new.cell.barcodes
  combined.matrix <- info %>% group_by(row) %>% dplyr::summarise(sum_aligned_reads_per_cell=sum(aligned_reads_per_cell, na.rm = TRUE), 
                                                          sum_aligned_umis_per_cell = sum(aligned_umis_per_cell, na.rm = TRUE))
  return(combined.matrix)
}

process.raw.rds <- function(data.file, sampleID, id, quality_score_cutoff=0.9,
                            default_cell_cutoff=200, dt.barcodes, n6.barcodes,
                            threshold=0.05){
  # merge DT and N6 primers
  counts <- merge.dtn6(data.file, sampleID, dt.barcodes, n6.barcodes)
  reads_umis_info <- merge.dtn6.aligned.reads.umis.per.cell(data.file, sampleID, dt.barcodes=dt.barcodes,
                                                            n6.barcodes=n6.barcodes)
  merged_aligned_reads_per_cell <- reads_umis_info$sum_aligned_reads_per_cell
  names(merged_aligned_reads_per_cell) <- reads_umis_info$row
  merged_umis_per_cell <- reads_umis_info$sum_aligned_umis_per_cell
  names(merged_umis_per_cell) <- reads_umis_info$row
  cells_df <- PrepareLqCellsData(counts, aligned.reads.per.cell=merged_aligned_reads_per_cell)
  
  cells_number_manual <- list(min=threshold*nrow(cells_df), max=nrow(cells_df))
  scores <- ScoreQualityData(merged_umis_per_cell, cells_df, cells_number_manual)
  filtered_cells <- names(scores)[which(scores > quality_score_cutoff)]
  filtered_counts <- counts[, colnames(counts) %in% filtered_cells]
  
  # remove cells with umi counts smaller than default cutoff value
  filtered_cells <- names(Matrix::colSums(filtered_counts)[Matrix::colSums(filtered_counts) > default_cell_cutoff])
  if(length(filtered_cells) == 1){
    filtered_counts <- as.matrix(filtered_counts[, colnames(filtered_counts) %in% filtered_cells])
    colnames(filtered_counts) <- filtered_cells
  } else{
    filtered_counts <- filtered_counts[, colnames(filtered_counts) %in% filtered_cells]
  }
  
  #push to seurat object
  seurat.obj <- CreateSeuratObject(counts = filtered_counts, project = id)
  # filter out > 200 < 7500 genes
  if(length(filtered_cells) > 1){
    seurat.obj <- subset(x = seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500)
  }
  
  original_counts <- counts
  save(seurat.obj,        
       original_counts, scores, merged_umis_per_cell, cells_df, cells_number_manual, 
       merged_aligned_reads_per_cell, reads_umis_info,
       file = paste0(id, ".seurat.filtered.RData")) 
  
  #return filtered count
  return(seurat.obj@assays$RNA@counts)
}

#############
# Parameters
args <- commandArgs(trailingOnly=TRUE)
# dropest rds file
data.file <- args[1]
# seurat object file name
id <- args[2]
# pagoda object file name
id2 <- args[3]

# sample table path
link.table <- read.table(args[4], sep="\t", header=TRUE)

# directory of R1_dTN6_pairs.txt file
paired.file <- args[5]

# QC filtering parameters
# remove cells with UMI_count smaller than 200
default_cell_cutoff <- 200
quality_score_cutoff <- 0.9
min_ncells <- 100

# processing id
#current project ID
sampleID <- gsub(".*/", "", id)
paired.barcodes <- read.table(paired.file,
                          header = F, stringsAsFactors = F)
dt.barcodes <- paired.barcodes[, 1]
n6.barcodes <- paired.barcodes[, 2]
merged_dtn6_cm_filtered <- process.raw.rds(data.file, sampleID, id, quality_score_cutoff=quality_score_cutoff,
                                           default_cell_cutoff=default_cell_cutoff,
                                           dt.barcodes=dt.barcodes, n6.barcodes=n6.barcodes)

# DT vs N6 
cm <- readRDS(data.file)$cm
cm_long_barcode <- colnames(cm)
cm_short8_barcode <- sapply(cm_long_barcode, function(c){
  str_sub(c,-8,-1)
})


cm_short8_barcode_inDT <- cm_short8_barcode %in% dt.barcodes
cm_short8_barcode_inN6 <- cm_short8_barcode %in% n6.barcodes

cm_withDT <- cm[, cm_short8_barcode_inDT]
cm_withN6 <- cm[, cm_short8_barcode_inN6]

DT_total_UMI_Number <- sum(cm_withDT)
N6_total_UMI_Number <- sum(cm_withN6)
DTvsN6 <- DT_total_UMI_Number/N6_total_UMI_Number

original_counts <- merge.dtn6(data.file, sampleID, dt.barcodes, n6.barcodes)
# save for pagoda
save(merged_dtn6_cm_filtered, original_counts, DT_total_UMI_Number,
        N6_total_UMI_Number, DTvsN6, file=paste0(id2, ".filtered.RData")) 

# generate p2 object - quick clustering
if(ncol(merged_dtn6_cm_filtered) > min_ncells){
  RNAp2 <- Pagoda2$new(merged_dtn6_cm_filtered, log.scale=TRUE, n.cores=5)
  RNAp2$adjustVariance(plot=F, gam.k=10)
  RNAp2$calculatePcaReduction(nPcs=30, n.odgenes=3e3)
  RNAp2$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine')
  RNAp2$getKnnClusters(method=infomap.community, type='PCA')
  RNAp2$getEmbedding(type='PCA', embeddingType='tSNE', perplexity=1, verbose=F)
  saveRDS(RNAp2, paste(id2, "p2.rds", sep="."))
}

