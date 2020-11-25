#' dimension reduction and visulization functions

#' caculate jaccard index matrix from an existing object
#' @param existing object: default pagoda2
#' @param max.cell max number of rows downsampled for calculating jaccard matrix
#' @param max.peaks max number of peaks, above it will do sampling
#' @export
RunJD <- function(obj, max.peaks=100000, max.cell=10000, seed.use=10){
  pmat <- obj[["pmat"]]
  if(is.null(pmat)){
    stop("there is no pmat slot in the object, please run addAtacObj and add peak matrix")
  }

  if((length(Matrix::rowMeans(pmat))) == 0L){
    stop("peak matrix is empty")
  }

  # if matrix is not binarized
  if((max(pmat)) > 1L){
    # binarized matrix
    #pmat[pmat > 0] <- 1
    pmat@x[pmat@x > 0] <- 1
  }
  p1 <- Matrix::rowMeans(pmat)

  if(ncol(pmat) > max.peaks){
    # downsampling peaks
    max.peaks <- min(max.peaks, ncol(pmat))
    col.covs <- log(Matrix::colSums(pmat)+1, 10);
    col.covs.dens <- density(col.covs, bw = 'nrd', adjust = 1)
    samplingProb <- 1 / (approx(col.covs.dens$x, col.covs.dens$y, xout=col.covs)$y + .Machine$double.eps)
    set.seed(seed.use)
    pid <- sort(sample(seq(col.covs), size=max.peaks, prob=samplingProb))
    pmat <- pmat[, pid]

    # remove the colums/rows have 0 coverage in the matrix
    pmat <- pmat[, which(log1p(Matrix::colSums(pmat)) > 0)]
  }

  if(max.cell < nrow(pmat)){
    max.cell <- min(max.cell, nrow(pmat))
    row.covs <- log(Matrix::rowSums(pmat)+1, 10)
    row.dens <- density(x = row.covs, bw = 'nrd', adjust = 1)
    samplingProb <- 1 / (approx(row.dens$x, row.dens$y, xout = row.covs)$y + .Machine$double.eps)
    set.seed(seed.use)
    pid <- sort(sample(x = seq(row.covs), size = max.cell, prob = samplingProb))
    pmat.ref = pmat[pid, ]
    p2 <- p1[pid]
  }else{
    pmat.ref <- pmat
    p2 <- p1
  }

  if(any(Matrix::rowSums(pmat) == 0)){
    print("removing zero rows in the matrix")
    pmat <- pmat[which(Matrix::rowSums(pmat)>0), ]
    pmat.ref <- pmat.ref[which(Matrix::rowSums(pmat.ref)>0), ]
  }

  jmat = calJaccard(pmat, pmat.ref)
  obj[["jmat"]]$jmat <- jmat
  obj[["jmat"]]$p1 <- p1
  obj[["jmat"]]$p2 <- p2

  return(obj)
}


#' run dimension reduction
#' @param obj p2 object
#' @param matrix matrix for dimension reduction (pmat, npmat, jmat, njmat) for ATAC data
dimReduct <- function(obj, matrix="npmat", method="svd", dm=50, scale.embedding=TRUE){

  dmatrix <- obj$counts
  if(matrix == "pmat"){
    dmatirx <- obj[["pmat"]]
  } else if(matrix == "npmat"){
    dmatrix <- obj[["npmat"]]
  } else if(matrix == "jmat"){
    dmatrix <- obj[["jmat"]]
  } else if(matrix=="njmat"){
    dmatrix <- obj[["njmat"]]
  }

  dm <- min(ncol(dmatrix), dm)
  embeddings <- list()

  if(method=="svd"){
    components <- irlba::irlba(A=dmatrix, nv=dm, work=dm+50)
    embeddings$loadings <- components$v
    embeddings$sdev <- components$d/sqrt(max(1, nrow(dmatrix)))

    embeddings$embed <- components$u
    if(scale.embedding){
      embeddings$embed <- t((t(embeddings$embed) - apply(embeddings$embed, 2, mean)) / apply(embeddings$embed, 2, sd))
    }
  }

  rownames(embeddings$embed) <- rownames(npmat)
  colnames(embeddings$embed) <- paste("svd", 1:ncol(embeddings$embed), sep="_")

  obj$reductions[["PCA"]] <- embeddings$embed
  #return(embeddings)
  return(obj)
}

#' plotComponentVariance based on pagoda obj
#' @param obj pagoda object
#' @param dm dimensions to plot
#' @export
plotComponentVariancePagoda <- function(obj, space='PCA', dm=20){
  if(is.null(obj$reductions[[space]])){
    stop("please run pca first and check if this is pagoda object")
  }

 reduction <- obj$reductions[[space]][, 1:dm]
 sd <- apply(reduction, 2, sd)
 data <- data.frame(PC=names(sd), sd=sd)
 data$PC <- factor(data$PC, levels = data$PC)

 ggplot2::ggplot(data=data, ggplot2::aes(x=PC, y=sd, group=1)) +
   ggplot2::geom_line() +
   ggplot2::geom_point(shape=16, size=2, alpha=0.8) +
   ggplot2::ylab('sd') + ggplot2::xlab('component number') +
   theme_bw()+
   theme(axis.text.x = element_text(angle = 90, hjust = 1))

}


#' run umap on PCA or other dm space
#' @param obj pagoda object
#' @param dims number of significant PCs for reduction (e.g. 1:20)
#' @export
quickUmap <- function(obj, space="PCA", dims=NULL, n.neighbors=30L, n.components=2L,
                      metric='cosine', n.epochs=NULL, learning.rate=1.0, min.dist=0.3,
                      spread=1.0, set.op.mix.ratio=1.0, local.connectivity=1L, repulsion.strength=1,
                      negative.sample.rate=5, a=NULL, b=NULL, uwot.sgd=FALSE, seed.use=22,
                      metric.kwds=NULL, angular.rp.forest=FALSE, verbose=TRUE){
  require(uwot)
  reduction <- obj$reductions[[space]]
  if(is.null(reduction)){
    stop("please run dimention reduction first")
  }

  if(!is.null(dims)){
    if(max(dims) > ncol(reduction)){
      warning("number of dims is larger than current variables, setting dims equal to it...")
      dims <- seq(min(dims), ncol(reduction), 1)
    }
    reduction <- reduction[, dims]
  }

  umap.embeddings <- uwot::umap(X = reduction,
    n_neighbors = as.integer(x = n.neighbors),
    n_components = as.integer(x = n.components),
    metric = metric,
    n_epochs = n.epochs,
    learning_rate = learning.rate,
    min_dist = min.dist,
    spread = spread,
    set_op_mix_ratio = set.op.mix.ratio,
    local_connectivity = local.connectivity,
    repulsion_strength = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    a = a, b = b,
    fast_sgd = uwot.sgd,
    verbose = verbose)
  rownames(umap.embeddings) <- rownames(reduction)

  obj$embeddings[[space]][["umap"]] <- umap.embeddings
  return(obj)
}

#' calculate high variable peaks based on normalized cell by peak matrix
#' @param obj pagoda2 object
#' @param percentileCutoff cutoff value for high variable peaks (0-1)
#' @export
hvgPeaks <- function(obj, percentileCutoff=NULL){
  npmat <- obj[["npmat"]]
  if(is.null(npmat)){
    stop("please normalize cell by peak matrix first")
  }

  peakCounts <- Matrix::colSums(npmat)
  # fit empirical cumulative distribution functon
  pcum <- ecdf(peakCounts)

  hvgpeaks <- data.frame(peak=colnames(npmat), peakCounts=peakCounts, percentile=pcum(peakCounts))
  hvgpeaks <- hvgpeaks[order(hvgpeaks$peakCounts, decreasing=TRUE), ]

  if(!is.null(percentileCutoff)){
    if(percentileCutoff > 1 | percentileCutoff < 0){
      stop("percentileCutoff is between 0 and 1")
    }
    hvgpeaks <- hvgpeaks[hvgpeaks$percentile>percentileCutoff, ]
  }
  return(hvgpeaks$peak)
}



#' plot features based on embedding for individual cluster
#' @param obj pagoda or related object
#' @export
plotFeature <- function(obj, clusterID=NULL, features="DAR", genename=NULL,
                        main=NULL, PvalCutoff=FALSE, cutoff=0.05){
  if(features == "DAR"){
    if(is.null(obj[["DAR"]])){
      stop("please run idenDAR first")
    }
    DAR <- obj[["DAR"]][[clusterID]]
    if(PvalCutoff){
      diffPeaks <- rownames(DAR[which(DAR$PValue < cutoff & DAR$logFC > 0), ])
    } else{
      diffPeaks <- rownames(DAR[which(DAR$adjPValue < cutoff & DAR$logFC > 0), ])
    }
    if(length(diffPeaks) < 3){
      return(NULL)
    }
    covs = Matrix::rowSums(obj[["pmat"]])
    vals = Matrix::rowSums(obj[["pmat"]][, which(colnames(obj[["pmat"]]) %in% diffPeaks)]) / covs
    vals.zscore = (vals - mean(vals)) / sd(vals)
    obj$plotEmbedding(type='PCA',embeddingType='tSNE',colors=vals.zscore,shuffle.colors=F,
                        mark.cluster.cex=1,alpha=0.7,main=main, quiet=T)
  } else if(features == "gene"){
    if(is.null(genename)){
      stop("please provide gene to visualize")
    }
    if(is.null(obj[["RNAcount"]])){
      stop("there is no correspondent RNA information in the object")
    }
    # matching
    RNAcount <- obj[["RNAcount"]]
    pmat <- obj[["pmat"]]

    RNAcount <- RNAcount[rownames(RNAcount) %in% rownames(pmat), ]
    if(is.null(main)){
      obj$plotEmbedding(type='PCA',embeddingType='tSNE', colors=RNAcount[, genename],
                        shuffle.colors=F,mark.cluster.cex=1,alpha=0.6, main=genename, quiet=T)
    } else{
      obj$plotEmbedding(type='PCA',embeddingType='tSNE', colors=RNAcount[, genename],
                        shuffle.colors=F,mark.cluster.cex=1,alpha=0.6, main=main, quiet=T)
    }
  }
}

