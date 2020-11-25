#' sample background sequences based on mathed GC content distribution from peaks
#' @param features sets of features (peaks) for comparison (Grange object)
#' @param genome genome sequences
sampleBackgroundSeq <- function(obj, features, genome="BSgenome.Hsapiens.NCBI.GRCh38", nsample=NULL){
  require(Biostrings)

  peaks <- obj[["peaks"]]
  if(is.null(peaks)){
    stop("please provide peak when run addAtacObj function")
  }

  if(is.null(nsample)){
    nsample <- length(features)
  }

  if(genome == "BSgenome.Hsapiens.NCBI.GRCh38"){
    require(BSgenome.Hsapiens.NCBI.GRCh38)
    genomeSeq <- BSgenome.Hsapiens.NCBI.GRCh38
    seqnames(genomeSeq) <- paste("chr", seqnames(genomeSeq), sep="")
    #species code
    spcode <- 9606

    # remove uncharaterized chromosomes
    chrlist <- paste("chr", seq(1:22), sep="")
    chrlist <- c(chrlist, "chrMT", "chrX", "chrY")
    peaks <- peaks[as.character(peaks@seqnames) %in% chrlist]
  }

  # remove uncharaterized chromosomes
  chrlist <- paste("chr", seq(1:22), sep="")
  chrlist <- c(chrlist, "chrMT", "chrX", "chrY")
  features <- features[as.character(features@seqnames) %in% chrlist]

  # get background regions
  backgroundRegion <- GenomicRanges::setdiff(peaks, features)

  peakSeq <- Biostrings::getSeq(genomeSeq, backgroundRegion)
  peakSeqCount <- Biostrings::alphabetFrequency(peakSeq)
  # GC distribution for all peaks
  GCbackground <- apply(peakSeqCount, 1, function(r){
    round((r[2] + r[3])/(r[1]+r[2]+r[3]+r[4]), 2)
  })

  featureSeq <- Biostrings::getSeq(genomeSeq, features)
  featureSeqCount <- Biostrings::alphabetFrequency(featureSeq)
  GCfeatures <- apply(featureSeqCount, 1, function(r){
    round((r[2] + r[3])/(r[1]+r[2]+r[3]+r[4]), 2)
  })

  # fit distribution
  featureWeights <- rep(0, length(backgroundRegion))
  densityEst <- density(GCfeatures, kernel = "gaussian", bw = 1)
  weights <- approx(densityEst$x, densityEst$y, xout=GCbackground,
                    yright = 0.00001,
                    yleft = 0.00001)$y
  featureWeights <- weights
  sampledRegion <- backgroundRegion[sample(seq(length(backgroundRegion)), size = nsample,
                          prob = featureWeights, replace=FALSE)]
  return(sampledRegion)
}

#' identify enriched motifs for a feature set (peaks) based on query regions
#' @param targetRegion genomic region sets for do enrichment analysis (e.g. chr17:44069893-44071328)
#' @param queryRegion genomic regions for all features (peaks)
enrichMotifs <- function(obj, targetRegion, sampleMethod="match.GC", nsample=10000){
  options("scipen"=0.0001, "digits"=4)
  motifMatrix <- obj[["motifMatrix"]]$matrix
  motifname <- obj[["motifMatrix"]]$motifname
  if(is.null(motifMatrix)){
    stop("please calculate motif matrix first using addMotifMatrix")
  }

  targetMotifs <- motifMatrix[rownames(motifMatrix) %in% targetRegion, ]

  # transform to Grange object
  targetChrs <- sapply(strsplit(targetRegion, ":"), `[`, 1)
  targetGenomicRegion <- sapply(strsplit(targetRegion, ":"), `[`, 2)
  targetRegionGrange <- GRanges(targetChrs, IRanges(as.numeric(sapply(strsplit(targetGenomicRegion, "-"), `[`, 1)),
                                                    as.numeric(sapply(strsplit(targetGenomicRegion, "-"), `[`, 2))))
  # sample regions
  backgroundSeq <- sampleBackgroundSeq(obj, targetRegionGrange, nsample=nsample)
  # transform back to characters
  backgroundSeq <- paste0(as.character(seqnames(backgroundSeq)), ":",
                          start(backgroundSeq), "-", end(backgroundSeq))
  motifMatrix.bg <- motifMatrix[backgroundSeq, ]


  targetMotifCounts <- Matrix::colSums(targetMotifs)
  bgMotifCounts <- Matrix::colSums(motifMatrix.bg)

  observed <- targetMotifCounts / length(targetRegion) * 100
  sampled <- bgMotifCounts / length(backgroundSeq) * 100
  fold.enrichment <- round(observed / sampled, 2)

  # hypergeometric test
  pVals <- unlist(lapply(seq(targetMotifCounts), function(r){
    phyper(targetMotifCounts[r]-1, bgMotifCounts[r],
           n=nrow(motifMatrix.bg) - bgMotifCounts[r],
           k=length(targetRegion), lower.tail=FALSE)
  }))

  enrichedMotifs <- data.frame(motif=colnames(targetMotifs),
                               observedCount=targetMotifCounts,
                               backgroundCount=bgMotifCounts,
                               observedPercent=observed,
                               backgroundPercent=sampled,
                               fold.enrichment=fold.enrichment,
                               pvalue = pVals, name=motifname)
  enrichedMotifs <- enrichedMotifs[order(enrichedMotifs$pvalue), ]
  return(enrichedMotifs)
}

#' identify overrepresentative motifs in DAR regions for each cluster
#' @param obj pagoda2 object
#' @param PvalueCutoff pvalue cutoff for significant DAR regions
#' @export
enrichMotifsDAR <- function(obj, sampleMethod="match.GC", PvalueCutoff=0.05, peakCutoff=500,
                            nsample=10000, genome="BSgenome.Hsapiens.NCBI.GRCh38", use.pvalues=FALSE){
  if(is.null(obj[["DAR"]])){
    stop("please run idenDAR to identify differatial accessible regions first")
  }
  dar <- obj[["DAR"]]
  clusters <- names(dar)

  enrichedMotifCluster <- parallel::mclapply(seq(length(clusters)), function(r){
    darRegion <- dar[[r]]
    if(use.pvalues){
      sigDARregion <- rownames(darRegion[darRegion$PValue < PvalueCutoff, ])
      sigDARregion <- sigDARregion[1:min(length(sigDARregion), peakCutoff)]
    } else{
      sigDARregion <- rownames(darRegion[darRegion$adjPValue<PvalueCutoff, ])
      sigDARregion <- sigDARregion[1:min(length(sigDARregion), peakCutoff)]
    }
    if(length(sigDARregion)>5){
      enrichedMotifCluter <- enrichMotifs(obj, sigDARregion, sampleMethod="match.GC", nsample=nsample)
    } else{
      enrichedMotifCluter <- NULL
    }
    enrichedMotifCluter
  }, mc.cores=5)

  names(enrichedMotifCluster) <- clusters
  obj[["enrichedMotifs"]] <- enrichedMotifCluster
  return(obj)
}

#' plot top motifs for each cluster
#' @export
plotTopMotifs <- function(obj, ntop=5, genome="human"){
  require(ggseqlogo)
  require(ggplot2)
  require(JASPAR2018)
  require(TFBSTools)

  enrichedMotifs <- obj[["enrichedMotifs"]]
  if(is.null(enrichedMotifs)){
    stop("please run enrichMotifsDAR first")
  }

  if(genome=="human"){
    pwms <- TFBSTools::getMatrixSet(JASPAR2018, opts = list(species = 9606, all_versions = FALSE))
  }

  motifPlot <- lapply(seq(as.numeric(names(enrichedMotifs))), function(r){
    clMotifs <- enrichedMotifs[[r]]
    if(!is.null(clMotifs) & !is.null(nrow(clMotifs[2]))){
      motifs <- rownames(clMotifs)[1:ntop]
      motifs.pwm <- pwms[names(pwms) %in% motifs]
      motif.profiles <- lapply(motifs.pwm, function(m){
        m@profileMatrix})
      motif.names <- unlist(lapply(motifs.pwm, function(m){
        m@name}))
      names(motif.profiles) <- motif.names
      ggseqlogo(motif.profiles, ncol=1) + ggplot2::ggtitle(paste("cluster", r, sep=" "))
    }
  })
  # remove clusters that have null element
  motifPlot <- plyr::compact(motifPlot)
  cowplot::plot_grid(plotlist=motifPlot, nrow=3)
}
