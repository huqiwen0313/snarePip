# sample-level QC feature generation
# usage: Rscript ATACpath sampleID RNAdir (NULL if not link)
# Qiwen Hu - 2020

library(SnapATAC)
library(pagoda2)
library(fastqcr)
library(snarePip)
library(GenomicRanges)
library(conos)
library(data.table)
library(DropletUtils)
library(dropestr)
suppressMessages(library(org.Hs.eg.db))

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
      fragment <- fragments[which(fragments[, 4] == as.character(barcodes[r])), ]
      fragmentUnique <- fragment[fragment$V5==1, ]
      fragment.region <- GRanges(as.character(fragment$V1), IRanges(as.numeric(fragment$V2), as.numeric(fragment$V3)))
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
    fragment.summary$peak.ratio <- round(fragment.summary$peak_overlap/fragment.summary$abundance, 3)
    fragment.summary$promoter.ratio <- round(fragment.summary$promoter_overlap/fragment.summary$uniqueFrag, 3)
    fragment.summary$blacklist.ratio <- round(fragment.summary$blacklist_overlap/fragment.summary$abundance, 3)
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

# arguments
args = commandArgs(trailingOnly=TRUE)
# path for ATAC processed data
ATACpath <- args[1]
# sample ID
sampleID <- args[2]
RNAdir <- args[3]
reference_genome <- args[4]
assay <- args[5]
link.table <- read.table(args[6], sep="\t", header=TRUE)
ref.dir <- args[7]
chromsize <- args[8]
pythonPath <- args[9]

# min cell size
min.cells.size <- 500
tssenrichmentCutoff <- 0.15
min.fragments <- 1000
# minimum number of dual cells after filtering
dual.ncells <- 500

tissue <- snarePip:::getTissue(link.table, sampleID)
# update ATAC/RNA paths based on file structure
ATACpath <- file.path(ATACpath, assay, tissue, "samples")
RNAdir <- file.path(RNAdir, assay, tissue, "samples")

# mt genes
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


## QC for alignment
bam.dir <- file.path(ATACpath, sampleID, "Sample_output", "bam")
bam.file <- paste(sampleID, "bam", sep=".")

# generate peak by cell matrix
snap.dir <- file.path(ATACpath, sampleID, "Sample_output", "snap")
snap.filename <- paste(sampleID, "snap", sep=".")
x.sp = createSnap(
  file=file.path(snap.dir, snap.filename),
  sample=sampleID,
  num.cores=10
)

x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=5)
x.sp = addPmatToSnap(x.sp)
median_number_duplicate_ratio <- round(median(1 - (x.sp@metaData$UQ+1)/(x.sp@metaData$PP+1)), 2)

# get peak by cell matrix and filtering
system(paste("mkdir", file.path(ATACpath, sampleID, "Sample_output", "pmats"), sep= " "),
       intern=TRUE)
peaks.dir <- file.path(ATACpath, sampleID, "Sample_output", "pmats")
# get valid chromosomes
chrlist <- paste0("chr", c(1:22, "X", "Y", "M"))
peaks <- x.sp@peak
ids <- which(as.character(peaks@seqnames) %in% chrlist)
peaks <- peaks[ids]

# filtering peak matrix
pmat <- x.sp@pmat[, ids]
peaks <- paste(peaks, sep=":")
colnames(pmat) <- peaks
pmat <- t(pagoda2::gene.vs.molecule.cell.filter(t(pmat), 
                                                min.cell.size=min.cells.size, plot=FALSE))

# summarize percell QCs
# fragment files
fragment.dir <- file.path(ATACpath, sampleID, "Sample_output", "fragments")
sample.fragment <- data.table::fread(file.path(fragment.dir, paste(sampleID, "fragments.bed", sep=".")), 
                                     showProgress = FALSE)

#sample.fragment[, 4] <- toupper(sample.fragment[, 4])
pos <- which(sample.fragment[, 2] > sample.fragment[, 3])
if(length(pos) > 0){
  sample.fragment <- sample.fragment[-1*pos, ]
}
sample.fragment <- sample.fragment[sample.fragment$V4 %in% rownames(pmat), ]

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
obj.dir <- file.path(ATACpath, sampleID, "Sample_output", "objects")
#obj.dir <- file.path(dir, "objects")
#atacPagoda <- Pagoda2$new(t(pmat), log.scale=TRUE)
atacPagoda <- readRDS(file.path(obj.dir, paste(sampleID, "p2.obj.rds", sep=".")))
atacPagoda <- addAtacObj(atacPagoda, pmat=pmat, fragments=sample.fragment, 
                         promoterRegion=promoters, TSS = tssRegion, blacklist=blacklist,
                         exons=exonRegion, transcripts=transcriptRegion, enhancers=enhancerRegion)

# calculate other statistics
atacPagoda$misc[["cellStat"]] <- NULL
stat <- cellStatistics(atacPagoda)

# tss enrichment
atacPagoda <- TSSenrichment(atacPagoda)
tss_enrichment <- round(median(atacPagoda[["TSSenrichment"]]), 2)

cellstat <- atacPagoda$misc[["cellStat"]]
tssenrichment <- atacPagoda[["TSSenrichment"]]
filteredCells <- names(tssenrichment[tssenrichment > tssenrichmentCutoff])
cellstatFiltered <- cellstat[(cellstat$barcode %in% filteredCells) & cellstat$abundance > min.fragments, ]
# save pmat
pmat <- pmat[rownames(pmat) %in% cellstatFiltered$barcode, ]
saveRDS(pmat, file.path(peaks.dir, paste(sampleID, ".pmat.filtered.rds", sep="")))
atacPagoda[["pmat"]] <- pmat

# save snap object
x.sp <- x.sp[which(x.sp@barcode %in% rownames(pmat)), ]
saveRDS(x.sp, file.path(snap.dir, paste(snap.filename, "rds", sep=".")))


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

overlapRatio <- 0
if(!is.null(atacPagoda[["RNAp2obj"]])){
  RNAp2 <- atacPagoda[["RNAp2obj"]]
  atacPagoda[["RNAcount"]] <- RNAp2$counts
  RNAbarcodes <- rownames(RNAp2$counts)
  atacbarcodes <- rownames(atacPagoda[["pmat"]])
  overlapRatio <- round(length(intersect(atacbarcodes, RNAbarcodes))/min(length(atacbarcodes), length(RNAbarcodes)), 2)
} else{
if(!is.null(RNAdir)){
  # link RNA with ATAC
  RNAfile <- file.path(RNAdir, paste(sampleID, "Sample_output", 
                                     "obj", 
                                     paste(sampleID, "sample_matrix.rds", sep="."), sep="/"))

  # read filtered rds
  RNAobj <- readRDS(RNAfile)
  
  # push to pagoda2 and run quick clustering
  RNAp2 <- Pagoda2$new(RNAobj, log.scale=TRUE, n.cores=2)
  RNAp2$adjustVariance(plot=F, gam.k=10)
  RNAp2$calculatePcaReduction(nPcs=30, n.odgenes=3e3)
  RNAp2$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine')
  RNAp2$getKnnClusters(method=infomap.community, type='PCA')
  RNAp2$getEmbedding(type='PCA', embeddingType='tSNE', perplexity=50, verbose=F)
  
  RNAbarcodes <- rownames(RNAp2$counts)
  atacbarcodes <- rownames(atacPagoda[["pmat"]])
  
  overlapRatio <- round(length(intersect(atacbarcodes, RNAbarcodes))/min(length(atacbarcodes), length(RNAbarcodes)), 2)
  atacPagoda[["RNAcount"]] <- RNAp2$counts
  atacPagoda[["RNAp2obj"]] <- RNAp2
 }
}

if(stat[1, 2]>dual.ncells){
  sample_quality = "Pass"
} else{
  sample_quality = "Not pass"
}
qc.stat <- c(median_reads_per_cell, median_fraction_read_in_mitochondrial, stat[2, 2],
                          stat[3, 2], stat[4, 2], stat[5, 2], stat[6, 2], 
                          tss_enrichment, median_number_duplicate_ratio, stat[1, 2], overlapRatio, sample_quality)
measurement <- c("median_reads_per_cell", "median_fraction_read_in_mitochondrial", "total_usable_fragments",
                 "median_total_usable_fragments_per_cell", "median_unique_fragments_per_cell",
                 "median_fraction_read_in_promoters", "median_fraction_read_in_peaks", "tss_enrichment",
                 "median_number_duplicate_ratio", "number_of_cells_after_filtering", "fraction_of_dual_omic_cells", "sample quality")
qc.stat <- data.frame(measurement=measurement, stat=qc.stat)

# save pagoda object and output
saveRDS(atacPagoda, file.path(obj.dir, paste(sampleID, "p2.obj.rds", sep=".")))
qc.dir <- file.path(ATACpath, sampleID, "Sample_output", "QCs")
system(paste("mkdir", qc.dir, sep= " "),
       intern=TRUE)
write.table(qc.stat, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t", 
            file.path(qc.dir, paste(sampleID, "qc.txt", sep=".")))

#############################################
# generate dual-omics QCs and objects
atacDual.dir <- file.path(ATACpath, sampleID, "Sample_output","dual_omics")
rnaDual.dir <- file.path(RNAdir, sampleID, "Sample_output", "dual_omics")
if(!file.exists(atacDual.dir)){
  dir.create(atacDual.dir)
}
if(!file.exists(rnaDual.dir)){
  dir.create(rnaDual.dir)
}

if(!is.null(RNAdir)){
  dual.cells <- intersect(atacbarcodes, RNAbarcodes)
  atacP2.dual <- atacPagoda
  dual.pmat <- atacPagoda[["pmat"]][rownames(atacPagoda[["pmat"]]) %in% dual.cells, ]
  dual.cellstat <- cellstatFiltered[cellstatFiltered$barcode %in% dual.cells, ]
  dual.tssenrichment <- atacPagoda[["TSSenrichment"]][names(atacPagoda[["TSSenrichment"]]) %in% dual.cells]
  atacP2.dual[["TSSenrichment"]] <- dual.tssenrichment
  atacP2.dual[["pmat"]] <- dual.pmat
  atacP2.dual[["cellStat"]] <- dual.cellstat
  
  ########### RNA part
  RNAcount <- readRDS(file.path(RNAdir, paste(sampleID, "Sample_output", 
                                              "obj", 
                                              paste(sampleID, "sample_matrix.rds", sep="."), sep="/")))
  RNAcount.dual <- RNAcount[, colnames(RNAcount) %in% dual.cells]
  RNAcount.dual.copy <- RNAcount.dual
  # remove mit genes
  RNAcount.dual <- RNAcount.dual[-which(rownames(RNAcount.dual) %in% mit_genes), ]
  # save count matrix
  saveRDS(RNAcount.dual, file.path(rnaDual.dir, paste(sampleID, "sample_matrix_dual.rds", sep=".")))
  
  # p2 object
  rnaP2.dual <- p2proc(RNAcount.dual, n.cores=20, min.cells.per.gene=5, min.cell.size=200, n.odgenes=2e3, 
                       get.largevis=FALSE, make.geneknn=FALSE, perplexity=50, log.scale = TRUE, 
                       batch = NULL, nPcs = 30, k = 30, get.tsne = TRUE)
  # doublet score
  scrublescore <- snarePip:::GetScrubletScores(RNAcount.dual, min.molecules.per.gene=10, method="scrublet",
                                    pythonPath=pythonPath)
  doublet.detect <- snarePip:::GetScrubletScores(RNAcount.dual, min.molecules.per.gene=10, method="doubletDetection",
                                      pythonPath=pythonPath)
  rnaP2.dual[["doublet"]] <- list(scrublet=scrublescore, doubletDetection=doublet.detect)
  saveRDS(rnaP2.dual, file.path(rnaDual.dir, paste(sampleID, "p2.dual.rds", sep=".")))
  
  # p2 app
  ids <- unlist(lapply(mget(colnames(rnaP2.dual$counts), org.Hs.egALIAS2EG, ifnotfound=NA), function(x) x[1]))
  rids <- names(ids); names(rids) <- ids
  go.env <- list2env(eapply(org.Hs.egGO2ALLEGS,function(x) as.character(na.omit(rids[x]))))
  rnaP2.dual$makeGeneKnnGraph(n.cores=5) 
  p2web <- snarePip:::GetPagodaWebApp(rnaP2.dual, rnaP2.dual$clusters$PCA$multilevel, additional.metadata=list(), 
                           verbose=T, go.sets=NULL, go.env=go.env, test.pathways=TRUE)
  # save object
  p2web$serializeToStaticFast(file.path(rnaDual.dir, paste(sampleID, "p2APP.dual.bin", sep=".")), verbose=T)
  
  #RNA QCs
  aligned_UMIs_per_cell <- Matrix::colSums(RNAcount.dual.copy)
  median_UMIs_per_cell <- median(aligned_UMIs_per_cell)
  percent_read_mapping_to_mitochrondrial<- round(GetGenesetFraction_overall(RNAcount.dual.copy, mit_genes), 4)
  median_fraction_read_in_mitochondrial <- round(median(GetGenesetFraction(RNAcount.dual.copy, mit_genes)), 4)
  qc.stat <- c(median_UMIs_per_cell, median_fraction_read_in_mitochondrial,
               percent_read_mapping_to_mitochrondrial, sample_quality,
               paste(obj.dir,"/",sampleID, ".p2App.bin",sep=""))
  
  measurement <- c("median_UMIs_per_cell", "median_fraction_read_in_mitochondrial",
                   "percent_read_mapping_to_mitochrondrial", 
                   "Sample_quality", "p2_webApp_address")
  rna.qc.stat <- data.frame(measurement=measurement, stat=qc.stat)
  write.table(rna.qc.stat, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t", 
              file.path(rnaDual.dir, paste(sampleID, "rna.dual.qc.txt", sep=".")))
  
  ############### ATAC
  # ATAC QCs
  atac.dual.stat <- c(length(unique(dual.cellstat$barcode)), median(dual.cellstat$abundance, na.rm=TRUE),
            median(dual.cellstat$uniqueFrag, na.rm=TRUE),
            median(dual.cellstat$promoter.ratio, na.rm=TRUE),
            median(dual.cellstat$peak.ratio, na.rm=TRUE))
  measurement <- c("total_dual_cells", "median_total_usable_fragments_per_cell",
                   "median_unique_fragments_per_cell", "median_fraction_read_in_promoters", "median_fraction_read_in_peaks")
  
  dual.stat <- data.frame(measurement=measurement, stat=atac.dual.stat)
  median_reads_per_cell <- median(dual.cellstat$abundance)
  median_fraction_read_in_mitochondrial <- median(dual.cellstat$chrM.ratio)
  atac.qc.stat.dual <- c(median_reads_per_cell, median_fraction_read_in_mitochondrial,
                         dual.stat[2, 2],  dual.stat[3, 2],  dual.stat[4, 2], 
                         dual.stat[5, 2], round(mean(dual.tssenrichment), 2), dual.stat[1, 2], sample_quality)
  measurement <- c("median_reads_per_cell", "median_fraction_read_in_mitochondrial",
                   "median_total_usable_fragments_per_cell", "median_unique_fragments_per_cell",
                   "median_fraction_read_in_promoters", "median_fraction_read_in_peaks", "mean_tss_enrichment",
                   "number of dual cells","sample quality")
  atac.qc.dual <- data.frame(measurement=measurement, stat=atac.qc.stat.dual)
  write.table(atac.qc.dual, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t", 
              file.path(atacDual.dir, paste(sampleID, "dual.qc.txt", sep=".")))
  
  # save pmat
  saveRDS(dual.pmat, file.path(atacDual.dir, paste(sampleID, "pmat_dual.rds", sep=".")))
  
  # p2 object
  atacP2.dual <- RunJD(atacP2.dual)
  atacP2.dual <- normJD(atacP2.dual)
  atacP2.dual$counts <- Matrix::Matrix(atacP2.dual[["jmat"]]$jmat, sparse=TRUE)
  atacP2.dual$adjustVariance(plot=F, gam.k=10)
  atacP2.dual$counts <- Matrix::Matrix(atacP2.dual[["njmat"]], sparse=TRUE)
  atacP2.dual$misc$odgenes <- colnames(atacP2.dual[["njmat"]])
  atacP2.dual$calculatePcaReduction(nPcs=50, n.odgenes=3e6)
  
  atacP2.dual$makeKnnGraph(k=15, type='PCA', center=T, distance='cosine')
  
  # differential accessibility regions
  if(nrow(dual.pmat) > 100){
    # get RNA clusters
    clusters <- rnaP2.dual$clusters$PCA$multilevel
    atacP2.dual <- idenDAR(atacP2.dual, clusters=as.factor(clusters), method="knn", 
                          test.method="lrt", bcv=0.4, n.cores=5)
  
  # motif analysis
    atacP2.dual <- addMotifMatrix(atacP2.dual, genome="BSgenome.Hsapiens.NCBI.GRCh38")
    atacP2.dual <- enrichMotifsDAR(atacP2.dual, sampleMethod="match.GC", PvalueCutoff=0.05, peakCutoff=500,
                                  nsample=10000, genome="BSgenome.Hsapiens.NCBI.GRCh38", use.pvalues=TRUE)
    
  # cicero activities
    #input_cds <- createInputCDS(pmat=dual.pmat)
    #gene.act <- ciceroGeneActivity(input_cds, genome="hg38", chromsize)
    #atacP2.dual[["cicero"]] <- gene.act
    saveRDS(atacP2.dual, file.path(atacDual.dir, paste(sampleID, "p2_dual.rds", sep=".")))
  }
}

