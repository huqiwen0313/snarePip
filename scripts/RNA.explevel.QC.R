# script for RNA experiment-level QCs
# Xin Wang - 2020

library(fastqcr)
library(dropestr)
library(Matrix)
require(stringr)
library(DropletUtils)
library(Seurat)

GetGenesetFraction_overall <- function (count.matrix, genes){
  umi.counts <- sort(Matrix::colSums(count.matrix), decreasing = T)
  presented.mit.genes <- intersect(genes, rownames(count.matrix))
  genes.frac <- sum(Matrix::colSums(count.matrix[presented.mit.genes, 
                                                 names(umi.counts)]))/sum(umi.counts)
  return(genes.frac)
}

# define args
args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("please specify reference_genome names", call.=FALSE)
}

sample <- args[1]
reference_genome <- args[2]
id = args[3]
mt_path <- args[4]

## QC for fastq files
qc.sample <- paste0("./tmp/","FastQC/",sample, "_R1_filtered_fastqc.zip")
qc <- qc_read(file.path(qc.sample))
quality_scores <- qc$per_sequence_quality_scores

# calculate fastq QC statistics
number_of_raw_reads <- as.numeric(qc$basic_statistics$Value[4])
quality_score <- round(sum(quality_scores[, 1]*quality_scores[, 2])/sum(quality_scores[, 2]), 4)
as.data.frame(quality_scores[, 1]>=30)[,1]->TnF_bt30
percent_q30_score <- round(sum(quality_scores[TnF_bt30, 2])/sum(quality_scores[, 2]), 4)

## QC for STAR alignment
bam.log.file <- paste0("./tmp/", "alignment/", sample, ".Log.final.out")
alignmentStat <- read.delim(bam.log.file)
percentage_uniquely_mapped_reads <- as.character(alignmentStat[8,2])
percent_uniquely_mapped_reads <- as.numeric(strsplit(percentage_uniquely_mapped_reads, "%")[[1]])/100

as.numeric(as.character(alignmentStat[22,2]))->multi_loci_mapped
as.numeric(as.character(alignmentStat[24,2]))->tooMany_loci_mapped
as.numeric(as.character(alignmentStat[7,2]))->unique_mapped_reads
percent_read_mapping_to_genome <- round((multi_loci_mapped+tooMany_loci_mapped+unique_mapped_reads)/number_of_raw_reads, 2)
read_length_aligned <- as.character(alignmentStat[5,2])

# Dropest
rdata.file <- paste0("./tmp/","pagoda_RData/", sample, ".filtered.RData")
load(rdata.file)
cm <- original_counts
umi_counts <- sort(Matrix::colSums(cm), decreasing=T)
mit_genes<- read.table(mt_path)$V1

total_cell_number <- ncol(cm)
cm_filtered <- merged_dtn6_cm_filtered
cell_number_after_filtering <- ncol(cm_filtered) 
cell_remain_rate_after_filtering <- signif(as.numeric(ncol(cm_filtered))/as.numeric(total_cell_number),3)


# stat before filtering 
percent_read_mapping_to_mitochrondrial <- round(GetGenesetFraction_overall(cm, mit_genes),4)
aligned_UMIs_per_cell <- Matrix::colSums(cm)
median_UMIs_per_cell <- median(aligned_UMIs_per_cell)
median_fraction_read_in_mitochondrial <- median(GetGenesetFraction(cm, mit_genes))

# stat After filtering
cm_filtered -> dcm_filtered
if(ncol(cm_filtered)>2){
  percent_read_mapping_to_mitochrondrial_after_cell_filtering <- round(GetGenesetFraction_overall(dcm_filtered, mit_genes),4)
  aligned_UMIs_per_cell_filtered <- apply(cm_filtered, 2, sum)
  median_UMIs_per_cell_after_cell_filtering <- median(aligned_UMIs_per_cell_filtered)
  median_fraction_read_in_mitochondrial_after_cell_filtering <- median(GetGenesetFraction(dcm_filtered, mit_genes))
  
  qc.stat <- c(read_length_aligned, reference_genome, number_of_raw_reads, quality_score,
               percent_q30_score, percent_read_mapping_to_genome, unique_mapped_reads,
               percent_uniquely_mapped_reads, percent_read_mapping_to_mitochrondrial,
               median_UMIs_per_cell, median_fraction_read_in_mitochondrial, total_cell_number,
               cell_number_after_filtering,cell_remain_rate_after_filtering, percent_read_mapping_to_mitochrondrial_after_cell_filtering,
               median_UMIs_per_cell_after_cell_filtering, median_fraction_read_in_mitochondrial_after_cell_filtering,
               DT_total_UMI_Number, N6_total_UMI_Number, DTvsN6)
  
  measurement <- c("read_length_aligned", "reference_genome", "number_of_raw_reads",
                   "quality_score", "percent_q30_score", "percent_read_mapping_to_genome",
                   "unique_mapped_reads", "percent_uniquely_mapped_reads", "percent_read_mapping_to_mitochrondrial",
                   "median_UMIs_per_cell", "median_fraction_read_in_mitochondrial", "total_cell_number",
                   "cell_number_after_filtering","cell_remain_rate_after_filtering", "percent_read_mapping_to_mitochrondrial_after_cell_filtering",
                   "median_UMIs_per_cell_after_cell_filtering", "median_fraction_read_in_mitochondrial_after_cell_filtering",
                   "DT_total_UMI_Number", "N6_total_UMI_Number", "DTvsN6")
  
  qc.stat <- data.frame(measurement=measurement, stat=qc.stat)
} else{
  qc.stat <- c(read_length_aligned, reference_genome, number_of_raw_reads, quality_score,
               percent_q30_score, percent_read_mapping_to_genome, unique_mapped_reads,
               percent_uniquely_mapped_reads, percent_read_mapping_to_mitochrondrial,
               median_UMIs_per_cell, median_fraction_read_in_mitochondrial, total_cell_number,
               cell_number_after_filtering,cell_remain_rate_after_filtering, 0, 0, 0, 
               DT_total_UMI_Number, N6_total_UMI_Number, DTvsN6)
  
  measurement <- c("read_length_aligned", "reference_genome", "number_of_raw_reads",
                   "quality_score", "percent_q30_score", "percent_read_mapping_to_genome",
                   "unique_mapped_reads", "percent_uniquely_mapped_reads", "percent_read_mapping_to_mitochrondrial",
                   "median_UMIs_per_cell", "median_fraction_read_in_mitochondrial", "total_cell_number",
                   "cell_number_after_filtering","cell_remain_rate_after_filtering", "percent_read_mapping_to_mitochrondrial_after_cell_filtering",
                   "median_UMIs_per_cell_after_cell_filtering", "median_fraction_read_in_mitochondrial_after_cell_filtering",
                   "DT_total_UMI_Number", "N6_total_UMI_Number", "DTvsN6")
  
  qc.stat <- data.frame(measurement=measurement, stat=qc.stat)
}

qc.dir <- paste0("./tmp/","QCs")
write.table(qc.stat, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t", 
            file.path(qc.dir, paste(sample, "qc.txt", sep=".")))
