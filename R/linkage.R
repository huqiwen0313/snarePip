#' asssign peak to its nearest gene
#' @param use.pvalue TRUE, use pvalue to identify DAR region - default FALSE
#' @param cutoff, adjPvalue (Pvalue) cutoff
#' @param test.asso TRUE perform statistical test the significance of identified peak-gene pairs
#' @export
assignPeakToGene <- function(obj, use.pvalues=FALSE, cutoff=0.05, test.asso=TRUE, n.cores=10){
  DAR <- obj[["DAR"]]
  if(is.null(DAR)){
    stop("please find differential accessibe regions first")
  }

  darRegions <- lapply(1:length(DAR), function(r){
    if(use.pvalues){
      regions <- rownames(DAR[[r]][DAR[[r]]$PValue < cutoff, ])
    } else{
      regions <- rownames(DAR[[r]][DAR[[r]]$adjPValue < cutoff, ])
    }
    if(length(regions) > 0){
      return(data.frame(region=regions, cluster=r))
    }
  }) %>% plyr::compact() %>% dplyr::bind_rows()
  darRegions <- unique(darRegions)
  # remove duplicated regions
  darRegions <- darRegions[!duplicated(darRegions$region), ]

  # transfer to IGrangers
  darChrs <- sapply(strsplit(darRegions$region, ":"), `[`, 1)
  darGenomicRegion <- sapply(strsplit(darRegions$region, ":"), `[`, 2)
  darRegionGrange <- GRanges(darChrs, IRanges(as.numeric(sapply(strsplit(darGenomicRegion, "-"), `[`, 1)),
                                              as.numeric(sapply(strsplit(darGenomicRegion, "-"), `[`, 2))))
  overlapRegions <- rep("", length(darRegionGrange))

  # get annotated genomic regions
  promoter <- obj[["promoterRegion"]]
  exons <- obj[["exons"]]
  transcripts <- obj[["transcripts"]]
  enhancer <- obj[["enhancers"]]

  # record overlapping regions
  overlapPromoter <- GenomicRanges::findOverlaps(darRegionGrange, promoter, select="first")
  overlapRegions[which(!is.na(overlapPromoter))] <- "promoter"
  overlapExons <- GenomicRanges::findOverlaps(darRegionGrange, exons, select="first")
  overlapRegions[which(!is.na(overlapExons))] <- "exon"

  overlapTranscripts <- GenomicRanges::findOverlaps(darRegionGrange, transcripts, select="first")
  transcriptRegion <- which(!is.na(overlapTranscripts))
  exonRegion <- which(!is.na(overlapExons))
  intronRegion <- setdiff(transcriptRegion, exonRegion)
  overlapRegions[intronRegion] <- "intron"
  overlapRegions[which(overlapRegions == "")] <- "intergenic"

  # find nearest genes
  hits <- nearest(darRegionGrange, transcripts, select="all")
  # choose the first transcript in the same gene
  hits <- data.frame(daridx=hits@from, geneidx=hits@to)
  hits <- hits[!duplicated(hits$daridx), ]
  assoGeneCoor <- transcripts[hits[, 2]]
  darRegionGrange <- darRegionGrange[hits[, 1]]
  overlapRegions <- overlapRegions[hits[, 1]]
  darRegions <- darRegions[hits[, 1], ]

  DARassGenes <- data.frame(chr=darRegionGrange@seqnames, peakStart=darRegionGrange@ranges@start,
                            peakEnd=darRegionGrange@ranges@width + darRegionGrange@ranges@start - 1,
                            transcriptStart=assoGeneCoor@ranges@start,
                            transcriptEnd=assoGeneCoor@ranges@width + assoGeneCoor@ranges@start - 1,
                            genename=names(assoGeneCoor),
                            distanceFromTss=darRegionGrange@ranges@width + darRegionGrange@ranges@start -
                              assoGeneCoor@ranges@start,
                            loc=overlapRegions, cluster=darRegions$cluster)

  if(test.asso){
    if(is.null(obj[["RNAp2obj"]])){
      stop("please run p2 on RNA first when set test.asso = TRUE")
    }

    pmat <- obj[["pmat"]]
    RNAcount <- obj[["RNAp2obj"]]$counts
    pmat <- pmat[rownames(pmat) %in% rownames(RNAcount), ]
    RNAcount <- RNAcount[rownames(RNAcount) %in% rownames(pmat), ]
    pmat <- pmat[match(rownames(RNAcount), rownames(pmat)), ]
    pmat@x[pmat@x > 1] <- 1

    models <- parallel::mclapply(1:nrow(DARassGenes), function(r){
      if(length(which(colnames(RNAcount) == DARassGenes[r, ]$genename)) > 0){
        if(sum(RNAcount[, DARassGenes[r, ]$genename]) > 0){
          peaks.id = paste(DARassGenes[r, ]$chr, paste(DARassGenes[r, ]$peakStart, DARassGenes[r, ]$peakEnd, sep="-"),
                           sep=":")
          model <- t(as.data.frame(summary(stats::glm(y ~ x, data=data.frame(y=pmat[, peaks.id], x=RNAcount[, DARassGenes[r, ]$genename]),
                                           family=binomial(link='logit')))[["coefficients"]]["x",]))
          rownames(model) <- peaks.id
          return(model)
        }
       }
    }, mc.cores=n.cores) %>% plyr::compact()

    models <- do.call(rbind, models)
    models[models[,"z value"] < 0, 4] <- 1
    rownames(DARassGenes) <- paste(DARassGenes$chr, paste(DARassGenes$peakStart, DARassGenes$peakEnd, sep="-"),
                                   sep=":")
    DARassGenes <- merge(DARassGenes, models, by="row.names")[-1]
    DARassGenes$logPval = -log10(DARassGenes$`Pr(>|z|)`)
  }

  obj[["DARassGenes"]] <- DARassGenes
  return(obj)
}

#' Get gene expression matrix of associated marker genes
#' @param count raw count matrix - rows are cells, columns are genes
#' @param groups cluster annotation
#' @param genes genelist
#' @param scale TRUE scale the data
GetMarkerExprMatrix <- function(count, groups, genes, scale=TRUE, logscale=TRUE){
  require(ComplexHeatmap)

  count.assGenes <- count[, colnames(count) %in% genes]
  genes <- genes[genes %in% colnames(count)]

  count.assGenes <- count[, match(genes, colnames(count))]
  groups <- groups[names(groups) %in% rownames(count)] %>% droplevels
  celltypeSize <- table(groups)
  # caculate avg expression for each celltype
  markerGeneMatrix <- lapply(1:length(names(celltypeSize)), function(r){
    cells <- names(groups[groups==names(celltypeSize)[r]])
    Matrix::colSums(count.assGenes[rownames(count.assGenes) %in% cells, ]) / celltypeSize[r]
  }) %>% dplyr::bind_cols() %>% as.data.frame

  rownames(markerGeneMatrix) <- make.unique(colnames(count.assGenes))
  colnames(markerGeneMatrix) <- names(celltypeSize)
  if(logscale){
    markerGeneMatrix <- log1p(markerGeneMatrix)
  }

  # scaling
  if(scale){
    markerGeneMatrix <- t(apply(markerGeneMatrix, 1, function(r){
      (r - mean(r))/sd(r)
    }))
  }
  return(markerGeneMatrix)
}


#' plot heatmap of associated genes
#' @export
plotAssoGenes <- function(obj, pal=colorRampPalette(c('dodgerblue1','grey95','indianred1'))(1024),
                          distanceCutoff=-1000, use.pvalue.cutoff=0.05,
                          RNAcluster=NULL){
  DARassGenes <- obj[["DARassGenes"]]
  DARassGenes <- DARassGenes[order(DARassGenes$cluster), ]
  if(is.null(DARassGenes)){
    stop("please run assignPeakToGene first")
  }

  # get RNA object
  p2RNA <- obj[["RNAp2obj"]]
  if(is.null(p2RNA)){
    stop("please provide RNA pagoda2 object")
  }

  # classify peak into different categories
  # within 1kb of TSS
  if(use.pvalue.cutoff > 0){
    DARassGenes <- DARassGenes[DARassGenes$`Pr(>|z|)`<use.pvalue.cutoff, ]
  }

  assGenesClose <- DARassGenes[DARassGenes$distanceFromTss>distanceCutoff & DARassGenes$distanceFromTss<0, ]
  # get genes from RNA
  RNAcount <- p2RNA$misc$rawCounts
  #RNAcount <- p2RNA$counts
  assGenesClose <- assGenesClose[assGenesClose$genename %in% colnames(RNAcount), ]
  RNAcountClose <- RNAcount[, colnames(RNAcount) %in% assGenesClose$genename]

  assGenesCloseMat <- GetMarkerExprMatrix(RNAcountClose, groups=RNAcluster , assGenesClose$genename, scale=TRUE)
  rannot <- data.frame(cluster=factor(assGenesClose$cluster, levels=unique(assGenesClose$cluster)))
  rownames(rannot) <- make.unique(as.character(assGenesClose$genename))

  ha <- ComplexHeatmap::HeatmapAnnotation(df=rannot, which = "row")
  h1 <- ComplexHeatmap::Heatmap(assGenesCloseMat,
                          cluster_rows=FALSE, cluster_columns=FALSE,
                          row_names_gp = grid::gpar(fontsize = 5), col=pal,
                          column_title="Expression of associated genes, close(<1kb)", name="expression") + ha
  # distal
  DARassGenes <- obj[["DARassGenes"]]
  DARassGenes <- DARassGenes[order(DARassGenes$cluster), ]
  if(use.pvalue.cutoff > 0){
    DARassGenes <- DARassGenes[DARassGenes$`Pr(>|z|)`<use.pvalue.cutoff, ]
  }
  assGenesFar <- DARassGenes[which(DARassGenes$distanceFromTss<distanceCutoff & DARassGenes$distanceFromTss<0), ]
  RNAcountFar <- RNAcount[, colnames(RNAcount) %in% assGenesFar$genename]
  assGenesFarMat <- GetMarkerExprMatrix(RNAcountFar, RNAcluster, assGenesFar$genename, scale=TRUE)
  rannot <- data.frame(cluster=factor(assGenesFar$cluster))
  rownames(rannot) <- make.unique(as.character(assGenesFar$genename))

  ha <- ComplexHeatmap::HeatmapAnnotation(df=rannot, which = "row")
  h2 <- ComplexHeatmap::Heatmap(assGenesFarMat,
                          cluster_rows=FALSE, cluster_columns=FALSE,
                          row_names_gp = grid::gpar(fontsize = 5), col=pal,
                          column_title="Expression of associated genes, distal(>1kb)", name="expression") + ha
  #others
  DARassGenes <- obj[["DARassGenes"]]
  DARassGenes <- DARassGenes[order(DARassGenes$cluster), ]
  if(use.pvalue.cutoff > 0){
    DARassGenes <- DARassGenes[DARassGenes$`Pr(>|z|)`<use.pvalue.cutoff, ]
  }
  assGenesOther <- DARassGenes[which(DARassGenes$distanceFromTss>0), ]
  RNAcountOther <- RNAcount[, colnames(RNAcount) %in% assGenesFar$genename]
  assGenesOtherMat <- GetMarkerExprMatrix(RNAcountOther, RNAcluster, assGenesOther$genename, scale=TRUE)
  rannot <- data.frame(cluster=factor(assGenesOther$cluster))
  rownames(rannot) <- make.unique(as.character(assGenesOther$genename))

  ha <- ComplexHeatmap::HeatmapAnnotation(df=rannot, which = "row")
  h3 <- ComplexHeatmap::Heatmap(assGenesOtherMat,
                          cluster_rows=FALSE, cluster_columns=FALSE,
                          row_names_gp = grid::gpar(fontsize = 5), col=pal,
                          column_title="Expression of associated genes, other(exons/introns)", name="expression") + ha
  return(list(h1, h2, h3))
}




