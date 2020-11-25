#' creating input CDS from custom peaks
#' @param p2 pagoda 2 object
#' @param pmat cell by peak matrix
#' @export
createInputCDS <- function(p2=NULL, pmat=NULL){
  require(cicero)
  require(annotables)
  if(is.null(p2)){
    if(is.null(pmat)){
      stop("please provide cell by peak matrix either in p2 object or by specifying pmat in the function")
    }
  } else{
    pmat <- p2[["pmat"]]
  }

  indata <- pmat
  indata@x[indata@x > 0] <- 1
  indata <- Matrix::t(indata)
  rownames(indata) <- gsub(":|-", "_", rownames(indata))

  cellinfo <- data.frame(cells=colnames(indata))
  rownames(cellinfo) <- colnames(indata)

  # format peak info
  chrs <- sapply(strsplit(rownames(indata), "_"), `[`, 1)
  peakStart <- sapply(strsplit(rownames(indata), "_"), `[`, 2)
  peakEnd <- sapply(strsplit(rownames(indata), "_"), `[`, 3)
  peakinfo <- data.frame(chr=chrs, bp1=peakStart, bp2=peakEnd)
  peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
  row.names(peakinfo) <- peakinfo$site_name

  # make CDS
  fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
  pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
  input_cds <-  suppressWarnings(newCellDataSet(indata,
                                                phenoData = pd,
                                                featureData = fd,
                                                expressionFamily=VGAM::binomialff(),
                                                lowerDetectionLimit=0))

  input_cds@expressionFamily@vfamily <- "binomialff"
  input_cds <- monocle::detectGenes(input_cds)
  input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
  return(input_cds)
}

#' calculate cicero gene activity score
#' @param input_cds processed cicero object
#' @param geneomeSize data.frame contains chromosome sizes
#' @export
ciceroGeneActivity <- function(input_cds, genome="hg38", genomeSize){
  input_cds %<>% estimateSizeFactors(.) %>%
    reduceDimension(., max_components = 2, num_dim=6,
                               reduction_method = 'tSNE', norm_method = "none")

  tsne_coords <- Matrix::t(reducedDimA(input_cds))
  rownames(tsne_coords) <- row.names(pData(input_cds))
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)

  # calculate co-accessibility
  conns <- run_cicero(cicero_cds, genomeSize)
  conns[is.na(conns$coaccess), ]$coaccess <- 0

  # prepare gene annotation
  if(genome=="hg38"){
    gene_annotation <- annotables::grch38
  }

  gene_annotation %<>% dplyr::select(chr, start, end, strand, symbol, biotype, description) %>% as.data.frame()
  gene_annotation$chr <- paste("chr", gene_annotation$chr, sep="")

  pos <- subset(gene_annotation, strand == "1")
  pos <- pos[order(pos$start),]
  pos <- pos[!duplicated(pos$symbol),]
  pos$end <- pos$start + 1

  neg <- subset(gene_annotation, strand == "-1")
  neg <- neg[order(neg$start, decreasing = TRUE),]
  neg <- neg[!duplicated(neg$symbol),]
  neg$start <- neg$end - 1

  gene_annotation_sub <- rbind(pos, neg)

  gene_annotation_sub <- gene_annotation_sub[,c(1:3, 5)]
  chr <- paste("chr", c(1:22, "X", "Y"), sep="")
  gene_annotation_sub <- gene_annotation_sub[gene_annotation_sub$chr %in% chr, ]
  names(gene_annotation_sub)[4] <- "gene"
  input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)
  # Generate gene activity scores
  # generate unnormalized gene activity matrix
  unnorm_ga <- build_gene_activity_matrix(input_cds, conns, coaccess_cutoff=-5)
  # remove any rows/columns with all zeroes
  unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, !Matrix::colSums(unnorm_ga) == 0]
  # make a list of num_genes_expressed
  num_genes <- pData(input_cds)$num_genes_expressed
  names(num_genes) <- row.names(pData(input_cds))
  # normalize
  cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
  return(list(conns=conns, cicero_gene_activities=cicero_gene_activities))
}

