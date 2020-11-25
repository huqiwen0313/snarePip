#' Caculate the number of fragment that overlapped with peak region and return summary statistics
#'   (median_total_usable_fragments_per_cell, median_unique_fragments_per_cell)
#' @param obj an class object (default pagoda2)
#' @param fitModel Whether to estimate fragments assuming fragment contamination rate is 0.02
#' @return plot and summary statistics
#' @export
fragmentOverlapPeaks <- function(obj, plot=T, fitModel=FALSE, n.cores=10){
  if(length(which(names(obj) %in% c("peaks", "fragments")))<2){
    stop("no ATAC data found, please creat object first")
  }

  peaks <- obj[["peaks"]]
  fragments <- obj[["fragments"]]
  barcodes <- rownames(obj$counts)
  #fragments.region <- GRanges(fragments[, 1], IRanges(fragments[, 2], fragments[, 3]))

  fragmentOverlap <- parallel::mclapply(1:length(barcodes), function(r){
     fragment <- fragments[fragments[, 4] %in% barcodes[r], ]
     fragmentUnique <- fragment[fragment[, 5]==1, ]
     fragment.region <- GRanges(fragment[, 1], IRanges(fragment[, 2], fragment[, 3]))
     overlap <- sum(countOverlaps(fragment.region, peaks))
     c(overlap, nrow(fragment), nrow(fragmentUnique))
  }, mc.cores=n.cores)

  fragment.summary <- data.frame(overlap=unlist(lapply(fragmentOverlap, `[[`, 1)),
                                 abundance=unlist(lapply(fragmentOverlap, `[[`, 2)),
                                 uniqueFrag=unlist(lapply(fragmentOverlap, `[[`, 3)))
  fragment.summary$barcode <- barcodes
  # add fragment per cell statistics into misc
  if(length(which(names(obj) %in% c("misc"))) == 0){
    obj[["misc"]][["fragmentStat"]] <- fragment.summary
  } else{
    obj$misc[["fragmentStat"]] <- fragment.summary
  }

  if(plot){
    fragment.summary <- fragment.summary[order(fragment.summary$abundance, decreasing=TRUE), ]
    fragment.summary$rank <- seq(1, nrow(fragment.summary), 1)

    y.max <- max(fragment.summary$overlap) * 1.5
    y.min <- min(fragment.summary$overlap)

    gg.base <- ggplot2::ggplot(fragment.summary, aes(x=rank, y=overlap))
    p1 <- gg.base + ggplot2::scale_x_log10(expand = c(0, 0)) +
      ggplot2::scale_y_log10(limits = c(y.min, y.max), expand = c(0, 0)) +
      stat_smooth(aes(x = seq(length(unique(rank)))),
                  se = F, method = "lm", formula = y ~ poly(x, 10), size = 0.5)
    p1 <- p1 + ggplot2::annotation_logticks(short=ggplot2::unit(2, 'pt'), mid=ggplot2::unit(3, 'pt'),
                                            long=ggplot2::unit(4, 'pt'))
    p1 <- p1+ylab("Fragment Overlapping peaks") + xlab("cell rank") + ggtitle("Cells") + theme_bw()


    p2 <- ggplot2::ggplot(fragment.summary, aes(x=fragment.summary$abundance)) +
          geom_histogram(aes(y=..density..), colour="black", fill="white")+
          geom_density(alpha=.2, fill="#FF6666") + theme_bw() +
          xlab("Fragment per barcode") + ggtitle("Fragment Distribution")
    cowplot::plot_grid(plotlist=list(p1, p2))
  }
}

#' plot fragment overlap peaks
#' @export
fragmentOverlapPeaksPlot <- function(obj, plot=T, fitModel=FALSE){
  if(length(which(names(obj) %in% c("peaks", "fragments")))<2){
    stop("no ATAC data found, please creat object first")
  }
  if(is.null(obj$misc[["cellStat"]])){
    stop("please run cellStatistics first")
  }

  fragment.summary <- obj$misc[["cellStat"]]
  fragment.summary <- fragment.summary[order(fragment.summary$abundance, decreasing=TRUE), ]
  fragment.summary$rank <- seq(1, nrow(fragment.summary), 1)

  y.max <- max(fragment.summary$peak_overlap) * 1.5
  y.min <- min(fragment.summary$peak_overlap+1)

  gg.base <- ggplot2::ggplot(fragment.summary, aes(x=rank, y=peak_overlap))
  p1 <- gg.base + ggplot2::scale_x_log10(expand = c(0, 0)) +
      ggplot2::scale_y_log10(limits = c(y.min, y.max), expand = c(0, 0)) +
      stat_smooth(aes(x = seq(length(unique(rank)))),
                  se = F, method = "lm", formula = y ~ poly(x, 10), size = 0.5)
  p1 <- p1 + ggplot2::annotation_logticks(short=ggplot2::unit(2, 'pt'), mid=ggplot2::unit(3, 'pt'),
                                            long=ggplot2::unit(4, 'pt'))
  p1 <- p1+ylab("Fragment Overlapping peaks") + xlab("cell rank") + ggtitle("Cells") + theme_bw()


  p2 <- ggplot2::ggplot(fragment.summary, aes(x=fragment.summary$abundance)) +
      geom_histogram(aes(y=..density..), colour="black", fill="white") +
      geom_density(alpha=.2, fill="#FF6666") + theme_bw() +
      xlab("Fragment per barcode") + ggtitle("Fragment Distribution")
  cowplot::plot_grid(plotlist=list(p1, p2))

}

#' Compute TSS enrichment score for each cell
#' @param promoters GRanges object contains promoters regions (-1000 and +1000 flanking regions of a TSS)
#' @export
TSSenrichment <- function(obj, regions=NULL, flankbp=100, sampleRegion=T, nsamples=5000, n.cores=20){
  if(is.null(regions)){
    regions <- obj[["TSS"]]
  }
  posMat <- intPosMat(obj, regions, sampleRegion=T, nsamples=5000, n.cores=n.cores)

  #average read depth in the 100(flankbp) bps at each of the end flanks
  flanking.ave <- rowMeans(as.matrix(posMat[, c(1:flankbp, (ncol(posMat)-flankbp+1):(ncol(posMat)))]))
  # fill up 0 values
  flanking.ave[flanking.ave==0] <- mean(flanking.ave)
  # normalization
  normMat <- posMat / flanking.ave

  TSSscore <- rowMeans(as.matrix(normMat[, 501:1500]))
  obj[["TSSenrichment"]] <- TSSscore

  obj$misc[["TSS.enrichment.matrix"]] <- normMat
  return(obj)
}

#' plot TSSenrichment
#' @export
TSSenrichmentPlot <- function(obj, cutoff=2.0, plotGroup=TRUE){
  enrichmentMat <- obj$misc[["TSS.enrichment.matrix"]]
  enrichmentScore <- obj[["TSSenrichment"]]

  # plot entire distribution
  meanEnrichment <- colMeans(as.matrix(enrichmentMat))
  pos <- seq(-1000, 1001, 1)

  enrichmentDistr <- data.frame(pos=pos, meanEnrichment=meanEnrichment)
  p <- ggplot2::ggplot(data=enrichmentDistr, aes(x=pos, y=meanEnrichment)) +
    geom_line(color="cornflowerblue") + theme_bw() + ylab("Relative Enrichment") +
    xlab("Relative Position (bp from TSS)")

  if(plotGroup){
    hqCells <- names(enrichmentScore[enrichmentScore>cutoff])
    lqCells <- names(enrichmentScore[enrichmentScore<=cutoff])
    hqMeanEnrichment <- colMeans(as.matrix(enrichmentMat[rownames(enrichmentMat) %in% hqCells, ]))
    lqMeanEnrichment <- colMeans(as.matrix(enrichmentMat[rownames(enrichmentMat) %in% lqCells, ]))
    enrichmentDistrGroup <- data.frame(pos=pos, meanEnrichment=hqMeanEnrichment, group="High") %>%
      rbind(data.frame(pos=pos, meanEnrichment=lqMeanEnrichment, group="Low"))

    # plotting
    p <- ggplot(data=enrichmentDistrGroup, aes(x=pos, y=meanEnrichment, group=group)) +
      geom_line(aes(color=group)) + theme_bw() + ylab("Relative Enrichment") +
      xlab("Relative Position (bp from TSS)")
    p <- p + facet_grid(.~ group)
  }
  return(p)
}

#' plotting fragment distribution
#'  Scatterplot showing the number of fragments/barcode and the percent of fragments overlapping peaks
#' @export
plotFragmentScatter <- function(obj){
  if(is.null(obj$misc[["cellStat"]])){
    stop("Please run cellStatistics first")
  }
  fragmentStat <- obj$misc[["cellStat"]]
  percentage <- round(fragmentStat$peak_overlap / fragmentStat$abundance, 3)
  fragmentStat$percentage <- percentage
  fragmentStat$abundance_log <- log1p(fragmentStat$abundance)

  ggplot2::ggplot(fragmentStat, aes(x=abundance_log, y=percentage)) + geom_point() + theme_bw() +
    xlab("Fragments per barcode, log scale") + ylab("Fraction fragment overlapping with peaks")
}

#' prepare data frame for plotting
#' @export
PrepareFragmentPlot <- function(obj=NULL, counts=NULL, breaks, estimated.number=F){
  if (estimated.number) {
    print("estimation function has not implemented yet")
  }

  if(is.null(counts)){
    counts <- obj$misc[["cellStat"]]$abundance
  }
  fragmentCount <- log10(sort(counts, decreasing=T))

  h <- hist(fragmentCount, breaks=breaks, plot=F)
  h$breaks <- h$breaks[1:(length(h$breaks)-1)]

  y.mults <- 10 ** h$breaks
  y.label <- '#Fragments * #barcode'

  plot.df <- data.frame(breaks=h$breaks, y=h$counts * y.mults)
  if(estimated.number) {
    print("estimation function has not implemented yet")
  }

  res <- list(plot.df=plot.df, y.label=y.label, fragmentCount=fragmentCount)
  return(res)
}

#' @export
PlotFragmentHist <- function(obj=NULL, counts=NULL, breaks=100, title=NULL, estimated.number=F, alpha=0.6,
                                show.legend=T){
  if(is.null(counts) & is.null(obj)){
    stop("please provide fragment counts or run fragmentOverlapPeaks first")
  }
  if(is.null(counts)){
    counts <- obj$misc[["cellStat"]]$abundance
  }

  df <- PrepareFragmentPlot(counts=counts, breaks=breaks, estimated.number=estimated.number)
  ylabel <- df$y.label
  df <- df$plot.df

  bar.width <- df$breaks[2] - df$breaks[1]
  bar_aes <- ggplot2::aes()

  if (show.legend) {
    gg.theme <- ggplot2::theme(legend.position=c(0.01, 0.99),
                               legend.justification=c(0, 1),
                               legend.background=ggplot2::element_rect(fill=ggplot2::alpha('white', 0.7)))
  } else {
    gg.theme <- ggplot2::theme(legend.position='none')
  }

  gg <- ggplot2::ggplot(df, ggplot2::aes(x=breaks, y=y)) +
    ggplot2::geom_bar(bar_aes, stat='identity', col=I('black'), alpha=alpha, size=0.05) +
    ggplot2::labs(x='log10(#fragments in cell)', y=ylabel) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    gg.theme + df$fill.scalse + theme_bw()

  return(gg)
}

#' calculate QC statistics per each cell
#'  blacklist ratio, promoter ratio, fragment stat
#' @export
cellStatistics <- function(obj, overallStat=TRUE, n.cores=20){
  require(GenomicRanges)
  if(length(which(names(obj) %in% c("fragments")))<1){
    stop("no fragment file found, please creat object with fragment info first")
  }
  if(is.null(obj$misc[["cellStat"]])){
    fragments <- obj[["fragments"]]
    fragments <- fragments[fragments$V4 %in% rownames(obj[["pmat"]]), ]

    promoter_overlap <- 0
    blacklist_overlap <- 0
    peak_overlap <- 0

    peaks <- obj[["peaks"]]
    promoters <- obj[["promoterRegion"]]
    blacklist <- obj[["blacklist"]]

    barcodes <- unique(fragments$V4)
    print(length(barcodes))
    cellstatistics <- parallel::mclapply(1:length(barcodes), function(r){
      fragment <- fragments[fragments[, 4] == barcodes[r], ]
      fragmentUnique <- fragment[fragment[, 5]==1, ]
      fragment.region <- GRanges(fragment[, 1], IRanges(fragment[, 2], fragment[, 3]))
      chrM.ratio <- round(nrow(fragment[grep("chrM", fragment$V1), ]) / nrow(fragment), 3)

      if(!is.na(peaks)){
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
    stat <- c(length(barcodes), length(fragments$V5), median(fragment.summary$abundance, na.rm=TRUE),
              median(fragment.summary$uniqueFrag, na.rm=TRUE),
              median(fragment.summary$promoter.ratio, na.rm=TRUE),
              median(fragment.summary$peak.ratio, na.rm=TRUE))
    measurement <- c("total_cells", "total_usable_fragments", "median_total_usable_fragments_per_cell",
                     "median_unique_fragments_per_cell", "median_fraction_read_in_promoters", "median_fraction_read_in_peaks")
    totalStat <- data.frame(measurement=measurement, stat=stat)
  }
  return(totalStat)
}

#' plot saturation
#' @export
plotSaturation <- function(obj, fragment=NULL){
  fragments <- obj[["fragments"]]
  counts <- data.frame(count=log1p(fragments$V5))
  p <- ggplot2::ggplot(counts, aes(x=count)) +
    geom_histogram(color="darkblue", fill="lightblue") +
    theme_bw() +
    xlab("log Read Count") + ylab("Frequency")
  return(p)
}

#' caculate alignment statistics
#' @export
alignmentStat<- function(bamfile, minimap=TRUE, read.tag.names=F){
  bam <- read.bam.tags(bamfile)
  total_mapped_reads <- length(bam$rname)
  if(minimap){
    unique_mapped_reads <- length(which(bam$tag$s2==0))
  } else{
    stop("other alignment have not implemented yet")
  }

  percent_uniquely_mapped_reads <- round(unique_mapped_reads/total_mapped_reads, 2)
  read_mapping_to_mitochrondrial <- length(which(bam$rname == "chrM"))
  percent_read_mapping_to_mitochrondrial <- round(read_mapping_to_mitochrondrial/total_mapped_reads, 2)
  read_length_aligned <- round(mean(bam$qwidth), 2)
  mappingStat <- c(total_mapped_reads, percent_uniquely_mapped_reads, percent_read_mapping_to_mitochrondrial,
                   read_length_aligned, unique_mapped_reads)
  names(mappingStat) <- c("total_mapped_reads", "percent_uniquely_mapped_reads", "percent_read_mapping_to_mitochrondrial",
                          "read_length_aligned", "unique_mapped_reads")
  return(mappingStat)
}

#' plot the comparison between RNA and ATAC barcode
#' @param obj pagoda2 object
#' @export
plotSnareBarcodeComparison <- function(obj){
  if(is.null(obj[["RNAcount"]])){
    stop("Please provide related RNA count to the object")
  }
  if(is.null(obj[["pmat"]])){
    stop("please call peak first")
  }
  cellSize <- sort(Matrix::rowSums(obj[["RNAcount"]]), decreasing=T)
  cellrank <- seq(cellSize)
  atacbarcodes <- rownames(obj[["pmat"]])
  RNAbarcodeSort <- names(cellSize)
  #######################################
  # temporarily removing prefix
  #atacbarcodes <- gsub(".*_", "", atacbarcodes)
  #RNAbarcodeSort <- gsub(".*_", "", RNAbarcodeSort)

  overlapRate <- unlist(lapply(cellrank, function(r){
    round(length(intersect(RNAbarcodeSort[1:r], atacbarcodes))/length(RNAbarcodeSort[1:r]),2)
  }))
  data <- data.frame(rank=cellrank, overlapRate=overlapRate)
  ggplot2::ggplot(data=data, aes(x=cellrank, y=overlapRate, group=1)) +
    theme_bw() + xlab("Cell rank in RNA") +
    ylab("Fraction barcodes (ATAC) overlapping with RNA") + stat_smooth(se=F)
}

#' plot distribution of DAR overlapping with annotated regions
#' @param obj p2 object
#' @param use.pvalues use pvalue to identify DAR region
#' @param cutoff p value (adj pvalue) cutoff
#' @export
plotDARoverlap <- function(obj, use.pvalues=FALSE, cutoff=0.05){
  DAR <- obj[["DAR"]]
  if(is.null(DAR)){
    stop("please find differential accessibe regions first")
  }

  darRegions <- unlist(lapply(1:length(DAR), function(r){
    if(use.pvalues){
      regions <- rownames(DAR[[r]][DAR[[r]]$PValue < cutoff, ])
      if(length(regions) == 0){
        regions <- rownames(DAR[[r]][DAR[[r]]$pval < cutoff, ])
      }
    } else{
      regions <- rownames(DAR[[r]][DAR[[r]]$adjPValue < cutoff, ])
      if(length(regions) == 0){
        regions <- rownames(DAR[[r]][DAR[[r]]$qval < cutoff, ])
      }
    }
    return(regions)
  }))

  # transfer to IGrangers
  darChrs <- sapply(strsplit(darRegions, ":"), `[`, 1)
  darGenomicRegion <- sapply(strsplit(darRegions, ":"), `[`, 2)
  darRegionGrange <- GRanges(darChrs, IRanges(as.numeric(sapply(strsplit(darGenomicRegion, "-"), `[`, 1)),
                                                    as.numeric(sapply(strsplit(darGenomicRegion, "-"), `[`, 2))))

  # get annotated genomic regions
  promoter <- obj[["promoterRegion"]]
  exons <- obj[["exons"]]
  transcripts <- obj[["transcripts"]]
  enhancer <- obj[["enhancers"]]

  # caculate overlap
  totalDAR <- length(darRegionGrange)
  if(!is.null(promoter)){
    overlapPromoter <- length(unique(GenomicRanges::findOverlaps(promoter, darRegionGrange, select="first")))
  } else{
    overlapPromoter <- 0
  }

  if(!is.null(exons)){
    overlapExons <- length(unique(GenomicRanges::findOverlaps(darRegionGrange, exons, select="first")))
  } else{
    overlapExons <- 0
  }

  if(!is.null(transcripts)){
    overlapTranscripts <- length(unique(GenomicRanges::findOverlaps(darRegionGrange, transcripts, select="first")))
  } else{
    overlapTranscripts <- 0
  }

  if(!is.null(enhancer)){
    overlapEnhancer <- length(unique(GenomicRanges::findOverlaps(darRegionGrange, enhancer, select="first")))
  } else{
    overlapEnhancer <- 0
  }

  overlapIntrons <- overlapTranscripts - overlapExons
  intergenic <- length(darRegionGrange) - (overlapPromoter+overlapExons+overlapIntrons+overlapEnhancer)
  if(intergenic<0){
    intergenic <- 0
  }
  total <- overlapPromoter+overlapExons+overlapIntrons+overlapEnhancer+intergenic
  proportion <- data.frame(regions=c("promoter", "exons", "Introns", "enhancer"),
                           overlapRatio=c(round(overlapPromoter/total, 2),
                                          round(overlapExons/total, 2),
                                          round(overlapIntrons/total, 2),
                                          round(overlapEnhancer/total, 2)))
                                          #round(intergenic/total, 2)))

  ggplot2::ggplot(data=proportion, aes(x=regions, y=overlapRatio)) +
    geom_bar(stat="identity", position=position_dodge(), fill="steelblue") +
    geom_text(aes(label=overlapRatio), vjust=1.6, color="white",
              position = position_dodge(0.9), size=3.5) + theme_bw() +
    ggplot2::ggtitle("DAR overlapped with annotated regions")
}

