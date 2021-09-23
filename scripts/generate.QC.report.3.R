# this script is used to generate SNARE ATAC QC and analysis report for each sample 
# Qiwen Hu - 2020
# Usage: Rscript genreate.QC.report.R dir(directory for sample) sampleName RNAsampledir (NULL if no paired RNA)
library(pagoda2)
library(SnapATAC)
library(fastqcr)
library(snarePip)
library(GenomicRanges)
library(conos)

cellStatistics <- function(obj, overallStat=TRUE, n.cores=20){
  require(GenomicRanges)
  if(length(which(names(obj) %in% c("fragments")))<1){
    stop("no fragment file found, please creat object with fragment info first")
  }
  if(is.null(obj$misc[["cellStat"]])){
    fragments <- obj[["fragments"]]
    if(is.null(fragments)){
      stop("please provide fragments file in p2 to caculate overlap statistics")
    }
    
    fragments <- fragments[fragments$V4 %in% rownames(obj[["pmat"]]), ]

    promoter_overlap <- 0
    blacklist_overlap <- 0
    peak_overlap <- 0

    peaks <- obj[["peaks"]]
    promoters <- obj[["promoterRegion"]]
    blacklist <- obj[["blacklist"]]

    barcodes <- unique(fragments$V4)
    cellstatistics <- parallel::mclapply(1:length(barcodes), function(r){
      fragment <- fragments[fragments[, 4] == barcodes[r], ]
      fragmentUnique <- fragment[fragment[, 5]==1, ]
      fragment.region <- GRanges(fragment[, 1], IRanges(fragment[, 2], fragment[, 3]))
      chrM.ratio <- round(nrow(fragment[grep("chrM", fragment$V1), ]) / nrow(fragment), 3)

      if(!is.na(peaks)){
        #peaks_overlap <- sum(countOverlaps(fragment.region, peaks))
        peaks_overlap <- length(which(GenomicRanges::countOverlaps(fragment.region, peaks, type="any")>0))
      }
      if(!is.na(promoters)){
        promoter_overlap <- length(which(GenomicRanges::countOverlaps(fragment.region, promoters)>0))
      }
      if(!is.na(blacklist)){
        blacklist_overlap <- length(which(GenomicRanges::countOverlaps(fragment.region, blacklist)>0))
      }
      c(peaks_overlap, promoter_overlap, blacklist_overlap, nrow(fragment), nrow(fragmentUnique), chrM.ratio)
    }, mc.cores=n.cores)

    fragment.summary <- data.frame(peak_overlap=unlist(lapply(cellstatistics, `[[`, 1)),
                                   promoter_overlap=unlist(lapply(cellstatistics, `[[`, 2)),
                                   blacklist_overlap=unlist(lapply(cellstatistics, `[[`, 3)),
                                   abundance=unlist(lapply(cellstatistics, `[[`, 4)),
                                   uniqueFrag=unlist(lapply(cellstatistics, `[[`, 5)),
                                   chrM.ratio=unlist(lapply(cellstatistics, `[[`, 6)))
    fragment.summary$peak.ratio <- round(fragment.summary$peak_overlap/fragment.summary$uniqueFrag, 3)
    fragment.summary$promoter.ratio <- round(fragment.summary$promoter_overlap/fragment.summary$uniqueFrag, 3)
    fragment.summary$blacklist.ratio <- round(fragment.summary$blacklist_overlap/fragment.summary$uniqueFrag, 3)
    fragment.summary$barcode <- barcodes
    obj$misc[["cellStat"]] <- fragment.summary
  } else{
    fragments <- obj[["fragments"]]
    fragments <- fragments[fragments$V4 %in% rownames(obj[["pmat"]]), ]
    barcodes <- unique(fragments$V4)
    fragment.summary <- obj$misc[["cellStat"]]
  }
  if(overallStat){
    stat <- c(length(barcodes), length(fragments$V5), median(fragment.summary$abundance, na.rm=TRUE), median(fragment.summary$uniqueFrag, na.rm=TRUE),
              median(fragment.summary$promoter.ratio, na.rm=TRUE), median(fragment.summary$peak.ratio, na.rm=TRUE))
    measurement <- c("total_cells", "total_usable_fragments", "median_total_usable_fragments_per_cell",
                     "median_unique_fragments_per_cell", "median_fraction_read_in_promoters", "median_fraction_read_in_peaks")
    totalStat <- data.frame(measurement=measurement, stat=stat)
  }
  return(totalStat)
}

# return RNA correspondent features
linkSamples <- function(linkTable, atacSampleName, assayfix=FALSE){
  if(assayfix){
    RNAtable <- linkTable[grep("SPL-R|SNARE2-R", linkTable$Experiment_ID), ]
    atacTable <- linkTable[grep("SPL-AC|SNARE2-AC", linkTable$Experiment_ID), ]
    # RNA table
    RNAtable$Experiment_ID <- gsub("[-|_]SPL-R", "", RNAtable$Experiment_ID)
    RNAtable$Experiment_ID <- gsub("[-|_]SNARE2-R", "", RNAtable$Experiment_ID)
    #RNAtable$Experiment_ID <- gsub("-SNARE2-R", "", RNAtable$Experiment_ID)
    
    # ATAC table
    atacTable$Experiment_ID <- gsub("-SPL-AC", "", atacTable$Experiment_ID)
    atacTable$Experiment_ID <- gsub("_SNARE2-AC", "", atacTable$Experiment_ID)
    atacTable$Experiment_ID <- gsub("-SNARE2-AC", "", atacTable$Experiment_ID)
    atacTable$Experiment_ID <- gsub("_P[0-9]*", "", atacTable$Experiment_ID)
  } else{
    RNAtable <- linkTable[linkTable$Type=="RNA", ]
    atacTable <- linkTable[linkTable$Type=="ATAC", ]
  }
  
  sampleinfo <- merge(RNAtable, atacTable, by=c("Library_ID", "Project"))
  atacSampleName <- gsub("_S.*", "", atacSampleName)
  atacSampleName <- gsub(".Ad1", "", atacSampleName)
  
  RNAsampleName <- strsplit(sampleinfo[sampleinfo$Experiment_ID.y == atacSampleName, 3], "_")[[1]]
  RNAsampleDate <- paste(RNAsampleName[1], RNAsampleName[2], sep="_")
  RNAsamplePair <- RNAsampleName[length(RNAsampleName)]
  return(c(RNAsampleDate, RNAsamplePair))
}

####
# define args and cutoffs
args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("please specify sample names", call.=FALSE)
}

dir <- args[1]
sample <- args[2]
RNAdir <- args[3]
reference_genome <- args[4]
assay <- args[5]
linktable <- read.table(args[6], sep="\t", header=TRUE)
ref.dir <- args[7]

# filtering cutoff
tssenrichmentCutoff <- 1.5
min.fragments <- 1000
min.umi <- 500

####
# start QC steps
## QC for fastq files
qc.dir <- file.path(dir, "fastqFiles")
qc.sample <- paste(sample, "fastqc.zip", sep="_")
qc <- qc_read(file.path(qc.dir, qc.sample))
quality_scores <- qc$per_sequence_quality_scores

# caculate fastq QC statistics
number_of_raw_reads <- as.numeric(qc$basic_statistics$Value[4])
quality_score <- round(sum(quality_scores[, 1]*quality_scores[, 2])/sum(quality_scores[, 2]), 2)
percent_q30_score <- round(sum(quality_scores[quality_scores[, 1]>=30, 2])/sum(quality_scores[, 2]), 2)

## QC for alignment
bam.dir <- file.path(dir, "bam")
bam.file <- paste(sample, "sorted.bam", sep=".")
alignmentStat <- alignmentStat(file.path(bam.dir, bam.file))
percent_uniquely_mapped_reads <- alignmentStat[2]
percent_read_mapping_to_mitochrondrial <- alignmentStat[3]
percent_read_mapping_to_genome <- round(alignmentStat[1]/number_of_raw_reads, 2)
read_length_aligned <- alignmentStat[4]
unique_mapped_reads <- alignmentStat[5]

## preparing QC for ATAC-seq data
# load and output peak matrix
sp.dir <- file.path(dir, "snap")
sp.sample <- paste(sample, "snap", sep=".")
x.sp = createSnap(
  file=file.path(sp.dir, sp.sample),
  sample=sample,
  num.cores=10
)
x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=5)
x.sp = addPmatToSnap(x.sp)
median_number_duplicate_ratio <- round(median(1 - (x.sp@metaData$UQ+1)/(x.sp@metaData$PP+1)), 2)

# save peak matrix
# get valid chromosomes
chrlist <- paste0("chr", c(1:22, "X", "Y", "M"))
peaks <- x.sp@peak
ids <- which(as.character(peaks@seqnames) %in% chrlist)
peaks <- peaks[ids]

# filtering peak matrix
pmat <- x.sp@pmat[, ids]

peaks.dir <- file.path(dir, "pmats")
peaks <- paste(peaks, sep=":")
#pmat <- x.sp@pmat
colnames(pmat) <- peaks
saveRDS(pmat, file.path(peaks.dir, paste(sample, "pmat.rds", sep=".")))

# brief filtering
pmat <- t(pagoda2::gene.vs.molecule.cell.filter(t(pmat), min.cell.size=min.umi, plot=FALSE))
# if too few cells, exit
if(nrow(pmat) < 5){
  qc.stat <- c(read_length_aligned, reference_genome, number_of_raw_reads, quality_score,
               percent_q30_score, percent_read_mapping_to_genome, unique_mapped_reads,
               percent_uniquely_mapped_reads, percent_read_mapping_to_mitochrondrial,
               "not Pass")
  measurement <- c("read_length_aligned", "reference_genome", "number_of_raw_reads",
                   "quality_score", "fraction_of_q30_score", "fraction_of_read_mapping_to_genome",
                   "unique_mapped_reads", "fraction_of_uniquely_mapped_reads", "fraction_of_read_mapping_to_mitochrondrial",
                   "sample quality")
  qc.stat <- data.frame(measurement=measurement, stat=qc.stat)
  qc.dir <- file.path(dir, "QCs")
  write.table(qc.stat, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t", 
              file.path(qc.dir, paste(sample, "qc.txt", sep=".")))
} else{
# summarize percell QCs
# fragment files
fragment.dir <- file.path(dir, "fragments")
sample.fragment <- read.table(file.path(fragment.dir, 
                                        paste(sample, "fragements.sort.bed", sep=".")))
sample.fragment[, 4] <- toupper(sample.fragment[, 4])
pos <- which(sample.fragment[, 2] > sample.fragment[, 3])
if(length(pos) > 0){
  sample.fragment <- sample.fragment[-1*pos, ]
}

# other statistics
# loading annotations
# blackregion
blacklist <- read.table(file.path(ref.dir, "hg38.blacklist.bed"), sep="\t")
blacklist <- GRanges(blacklist[, 1], IRanges(blacklist[, 2], blacklist[, 3]))

# promotors
promoters <- read.table(file.path(ref.dir, "hg38.promoters.bed"), sep="\t")
promoters <- GRanges(promoters[,1], IRanges(promoters[, 2], promoters[, 3]))

# read tss sites
tssSites <- read.table(file.path(ref.dir, "hg38.tss.bed"))
tssSites <- tssSites[tssSites[, 7]=="protein_coding", ]
flanking <- 1000 
tssSites[, 2] <- tssSites[, 2] - flanking
tssSites[, 3] <- tssSites[, 3] + flanking
tssRegion <- GRanges(tssSites[, 1], IRanges(tssSites[, 2], tssSites[, 3]))

# exon annotation
gene.annot <- read.table(file.path(ref.dir, "hg38.genes.gtf"), sep="\t")
gene.names <- vapply(strsplit(as.character(gene.annot[, 9]), "; "), `[`, 4, FUN.VALUE=character(1))
gene.names <- gsub("gene_name ", "", gene.names)
gene.annot$geneName <- gene.names
exons <- gene.annot[gene.annot[, 3] == "exon", ]
exonRegion <- GRanges(exons[, 1], IRanges(exons[, 4], exons[, 5]))
names(exonRegion) <- exons$geneName

# transcript annotation
transcripts <- gene.annot[gene.annot[, 3]=="transcript", ]
transcript.names <- vapply(strsplit(as.character(transcripts[, 9]), "; "), `[`, 4, FUN.VALUE=character(1))
transcriptRegion <- GRanges(transcripts[, 1], IRanges(transcripts[, 4], transcripts[, 5]))
names(transcriptRegion) <- transcripts$geneName

# enhancers
enhancers <- read.table(file.path(ref.dir, "hg38.enhancer.bed"), sep="\t")
enhancerRegion <- GRanges(enhancers[, 1], IRanges(enhancers[, 2], enhancers[, 3]))

# add ATAC-features into existed object
obj.dir <- file.path(dir, "objects")
atacPagoda <- Pagoda2$new(t(pmat), log.scale=TRUE)
atacPagoda <- addAtacObj(atacPagoda, pmat=pmat, fragments=sample.fragment, 
                         promoterRegion=promoters, TSS = tssRegion, blacklist=blacklist,
                         exons=exonRegion, transcripts=transcriptRegion, enhancers=enhancerRegion)

# calculate other statistics
stat <- cellStatistics(atacPagoda)

# tss enrichment
atacPagoda <- TSSenrichment(atacPagoda)
tss_enrichment <- round(median(atacPagoda[["TSSenrichment"]]), 2)

# briefly filter low quality cells (promotor ratio > 0.15)
#promotorRationCutoff <- 0.15
cellstat <- atacPagoda$misc[["cellStat"]]
tssenrichment <- atacPagoda[["TSSenrichment"]]
filteredCells <- names(tssenrichment[tssenrichment > tssenrichmentCutoff])
cellstatFiltered <- cellstat[(cellstat$barcode %in% filteredCells) & cellstat$abundance > min.fragments, ]
# get filtered pmat
pmat <- pmat[rownames(pmat) %in% cellstatFiltered$barcode, ]
saveRDS(pmat, file.path(peaks.dir, paste(sample, "pmat.filtered.rds", sep=".")))
atacPagoda[["pmat"]] <- pmat

stat <- c(length(unique(cellstatFiltered$barcode)), stat[2, 2], median(cellstatFiltered$abundance, na.rm=TRUE),
          median(cellstatFiltered$uniqueFrag, na.rm=TRUE),
          median(cellstatFiltered$promoter.ratio, na.rm=TRUE),
          median(cellstatFiltered$peak.ratio, na.rm=TRUE))
measurement <- c("total_cells", "total_usable_fragments", "median_total_usable_fragments_per_cell",
                 "median_unique_fragments_per_cell", "median_fraction_read_in_promoters", "median_fraction_read_in_peaks")
stat <- data.frame(measurement=measurement, stat=stat)
median_reads_per_cell <- median(cellstatFiltered$abundance)
median_fraction_read_in_mitochondrial <- median(cellstatFiltered$chrM.ratio)

# output txt format of QC report
# start generate analysis report object
if(nrow(pmat) > 100){
  atacPagoda <- RunJD(atacPagoda)
  atacPagoda <- normJD(atacPagoda)
  atacPagoda$counts <- Matrix::Matrix(atacPagoda[["jmat"]]$jmat, sparse=TRUE)
  atacPagoda$adjustVariance(plot=F, gam.k=10)
  atacPagoda$counts <- Matrix::Matrix(atacPagoda[["njmat"]], sparse=TRUE)
  atacPagoda$misc$odgenes <- colnames(atacPagoda[["njmat"]])
  atacPagoda$calculatePcaReduction(nPcs=50, n.odgenes=3e6)
  
  atacPagoda$makeKnnGraph(k=15, type='PCA', center=T, distance='cosine')
  p2graph <- atacPagoda$graphs$PCA
  clusters <- igraph::membership(conos::leiden.community(p2graph, resolution=2))
}

# differential accessibility regions
if(nrow(pmat) > 100){
  atacPagoda <- idenDAR(atacPagoda, clusters=as.factor(clusters), method="knn", 
                        test.method="lrt", bcv=0.4, n.cores=5)
}

# motif analysis
if(nrow(pmat) > 100){
  atacPagoda <- addMotifMatrix(atacPagoda, genome="BSgenome.Hsapiens.NCBI.GRCh38")
  atacPagoda <- enrichMotifsDAR(atacPagoda, sampleMethod="match.GC", PvalueCutoff=0.05, peakCutoff=500,
                            nsample=10000, genome="BSgenome.Hsapiens.NCBI.GRCh38", use.pvalues=TRUE)
}

# link RNA
overlapRatio <- 0
if(!(RNAdir == "NULL")){
  # link RNA with ATAC
  RNAsampleinfo <- linkSamples(linktable, sample)
  libID <- snarePip:::getLibID(linktable, sample)
  experiment_id <- RNAsampleinfo[1]
  tissue <- snarePip:::getTissue(linktable, sampleName=experiment_id)
  RNAfile <- file.path(RNAdir, paste(assay, tissue, "samples", experiment_id, "Experiment_output",
                                     "pagoda_RData",
                                     paste(libID, "filtered.RData", sep="."), sep="/"))
  # read filtered rds
  t <- load(RNAfile)
  RNAobj <- get(t[1])

  # push to pagoda2 and run quick clustering
  if(ncol(RNAobj) > 100){
    RNAp2 <- Pagoda2$new(RNAobj, log.scale=TRUE, n.cores=2)
    RNAp2$adjustVariance(plot=F, gam.k=10)
    RNAp2$calculatePcaReduction(nPcs=50, n.odgenes=3e3)
    RNAp2$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine')
    RNAp2$getKnnClusters(method=infomap.community, type='PCA')
    RNAp2$getEmbedding(type='PCA', embeddingType='tSNE', perplexity=10, verbose=F)
    
    RNAbarcodes <- rownames(RNAp2$counts)
    atacbarcodes <- rownames(atacPagoda[["pmat"]])
    
    overlapRatio <- round(length(intersect(atacbarcodes, RNAbarcodes))/min(length(atacbarcodes), length(RNAbarcodes)), 2)
    atacPagoda[["RNAcount"]] <- RNAp2$counts
    atacPagoda[["RNAp2obj"]] <- RNAp2
  } else{
    overlapRatio <- 0
  }
}

qc.stat <- c(read_length_aligned, reference_genome, number_of_raw_reads, quality_score,
             percent_q30_score, percent_read_mapping_to_genome, unique_mapped_reads,
             percent_uniquely_mapped_reads, percent_read_mapping_to_mitochrondrial,
             median_reads_per_cell, median_fraction_read_in_mitochondrial, stat[2, 2],
             stat[3, 2], stat[4, 2], stat[5, 2], stat[6, 2], 
             tss_enrichment, median_number_duplicate_ratio, stat[1, 2], overlapRatio)
measurement <- c("read_length_aligned", "reference_genome", "number_of_raw_reads",
                 "quality_score", "fraction_of_q30_score", "fraction_of_read_mapping_to_genome",
                 "unique_mapped_reads", "fraction_of_uniquely_mapped_reads", "fraction_of_read_mapping_to_mitochrondrial",
                 "median_reads_per_cell", "median_fraction_read_in_mitochondrial", "total_usable_fragments",
                 "median_total_usable_fragments_per_cell", "median_unique_fragments_per_cell",
                 "median_fraction_read_in_promoters", "median_fraction_read_in_peaks", "tss_enrichment",
                 "median_number_duplicate_ratio", "number_of_cells_after_filtering", "fraction_of_dual_omic_cells")
qc.stat <- data.frame(measurement=measurement, stat=qc.stat)

# save pagoda object and output
saveRDS(atacPagoda, file.path(obj.dir, paste(sample, "p2.obj.rds", sep=".")))
qc.dir <- file.path(dir, "QCs")
write.table(qc.stat, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t", 
            file.path(qc.dir, paste(sample, "qc.txt", sep=".")))
}
