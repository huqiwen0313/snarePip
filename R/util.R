#' read bam files and extract tags from minimap2
#' @param filename bam file
#' @import Rsamtools
#' @export
read.bam.tags <- function(filename, read.tag.names=F, fix.chromosome.names=F){
  if(!is.element("Rsamtools", installed.packages()[, 1])) {
    stop("Rsamtools Bioconductor package is now required for BAM file support. Please install")
  }

  ww <- c("flag","rname","pos","isize","strand","mapq","qwidth"); if(read.tag.names) { ww <- c(ww,"qname") };
  bam <- Rsamtools::scanBam(filename,param=Rsamtools::ScanBamParam(what=ww,tag="s2",
                                                                   flag=Rsamtools::scanBamFlag(isUnmappedQuery=FALSE)))[[1]]
  return(bam)
}


#' Generate cell position matrix in a given region based on middle position of a fragment
#' @param obj an R object (default:pagoda2)
#' @param region regions of interest (GRanges object)
#' @param cells which cell to include. Default: all cells in the object
#' @export
posMatrix <- function(obj, region, cells=NULL, fragments=NULL){
  if(is.null(obj[["fragments"]]) & is.null(fragments)){
    stop("please provide fragments information")
  }
  if(is.null(fragments)){
    fragments <- obj[["fragments"]]
  }
  if(is.null(cells)){
    cells <- rownames(obj[["pmat"]])
  }
  fragments <- fragments[fragments$V4 %in% cells, ]
  # extract cell fragments fall into a given region (based on middle position of fragment)
  fragments <- fragments %>% subset(V1==as.character(region@seqnames))
  cellFragments <- data.frame(pos=round(fragments[, 2]+abs(fragments[, 3]-fragments[, 2]+1)/2, 0)-
                                      start(region)+1,
                    cells=c(fragments[, 4]),
                    stringsAsFactors = FALSE) %>% subset(pos>0 & pos<=width(region))

  cellindex <- setNames(seq_along(cells), cells)

  # creat position matrix
  if(nrow(cellFragments) == 0){
    posMat <- sparseMatrix(
      i = NULL,
      j = NULL,
      dims = c(length(cells), width(region))
    )
  } else{
    posMat <- sparseMatrix(
      i = cellindex[cellFragments$cells],
      j = cellFragments$pos,
      x = 1,
      dims = c(length(cells), width(region))
    )
  }
  rownames(x = posMat) <- cells
  colnames(x = posMat) <- start(region):end(region)

  return(posMat)
}

#' generate integrated cell position matrix given sets of regions
#' @param regions regions(GRanges object) that used to generate cell position matrix
#' @return position matrix
#' @export
intPosMat <- function(obj, regions, cells=NULL, sampleRegion=T, nsamples=5000, n.cores=10){
  if(sampleRegion){
    regions <- sample(regions, nsamples)
  }

  posMat <- parallel::mclapply(seq_along(regions), function(r){
    posMatrix(obj, regions[r, ])
  }, mc.cores=n.cores) %>% purrr::reduce(`+`)
  return(posMat)
}

#' calculate jaccard index between two binary matrices
#' @param m1 matrix 1
#' @param m2 matrix 2
#' @export
calJaccard <- function(m1, m2){
  if(any(m1) > 1 | any(m2) > 1){
    stop("matrix contain values larger than 1, please binarize matrix")
  }

  m1.count.all <- Matrix::rowSums(m1)
  m2.count.all <- Matrix::rowSums(m2)

  if(length(which(m1.count.all==0))>0 | length(which(m2.count.all==0))>0){
    stop("please remove zero rows in the matrix")
  }

  crosspd <- Matrix::tcrossprod(m1, m2)
  #crosspd <- fast_multiply(m1, m2)
  m1i <- replicate(ncol(crosspd), m1.count.all)
  m2i <- replicate(nrow(crosspd), m2.count.all)

  jmat = as.matrix(crosspd / ((m1i + t(m2i)) - crosspd))
  #jmat = crosspd / ((m1i + t(m2i)) - crosspd)
  rownames(jmat) <- rownames(m1)
  colnames(jmat) <- rownames(m2)
  return(jmat)
}

#' normalize peak matrix based on TF-IDF and other algorithms
#' @param obj p2 object
#' @param method TfIDF: LSI TF-IDF normalization in Stuart & Butler et al. 2019
#'         log-TfIDF: The log-TF method
#'         TF: IDF normalization no TF
#' @return normalized peak matrix
#' @export
normPmatrix <- function(obj, method="TfIDF", scale.factor=1e4){
  pmat <- t(obj[["pmat"]])

  if(is.null(pmat)){
    stop("please provide peak by cell matrix into pagoda object")
  }

  message("normalize peak by cell matrix")
  npeaks <- Matrix::colSums(pmat)
  tf <- Matrix::tcrossprod(pmat, Matrix::Diagonal(x=1/npeaks))

  if(method == "TfIDF"){
    idf <- ncol(pmat) / Matrix::rowSums(pmat)
  } else if(method == "log-TfIDF"){
    idf <- log((1+(ncol(pmat)/rowSums(pmat))), base=2)
  } else if(method == "IDF"){
    tf <- pmat
  }

  npmat <- Diagonal(n=length(idf), x=idf) %*% tf

  if (method=="TfIDF"){
    npmat <- log1p(npmat * scale.factor)
  }

  obj[["npmat"]] <- t(npmat)
  return(obj)
}

#' normalize jaccard index matrix based on fang et al
#'   doi: https://doi.org/10.1101/615179
#' @param obj default: pagoda object
#' @return obj with normalized jaccard matrix
#' @export
normJD <- function(obj, method="residual", center=TRUE){
  jmat <- obj[["jmat"]]$jmat
  #pmat <- obj[["pmat"]]

  #if((max(pmat)) > 1L){
    # binarized matrix
    #pmat[pmat > 0] <- 1
  #  pmat@x[pmat@x > 0] <- 1
  #}

  if(is.null(jmat)){
    stop("please calculate jaccard index matrix first")
  }

  p1 = obj[["jmat"]]$p1
  p2 = obj[["jmat"]]$p2

  # fit polynormial regression model
  # expected value based on read depth
  cspd <- tcrossprod(p1, p2)
  dsum <- matrix(rep(p1, each=length(p2)), ncol=length(p2), byrow=TRUE) +
      matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
  expJaccard <- cspd/(dsum - cspd)
  scale.factor <- mean(jmat / expJaccard)

  diag(jmat) <- scale.factor * diag(expJaccard)

  data <- data.frame(x=c(expJaccard), y=c(jmat))

  # fit model
  model <- lm(y ~ x + I(x^2), data)
  preds = predict(model, data.frame(x=c(expJaccard)), se.fit = TRUE)

  if(method == "zscore"){
    nmatrix = (c(jmat) - preds$fit) / (preds$se.fit);
  }else if(method == "residual"){
    nmatrix = c(jmat) -  preds$fit;
  }
  nmat = matrix(nmatrix, nrow(expJaccard), ncol(expJaccard))
  rownames(nmat) <- rownames(jmat)
  colnames(nmat) <- colnames(jmat)

  if(center){
    nmat = t(scale(t(nmat), center=TRUE, scale=TRUE))
  }

  obj[["njmat"]] <- nmat
  return(obj)
}

#' sampling background cells
#' @param obj pagoda2 object
#' @param posCell position of positive cells
#' @param method, sample cells based on nearest cells or random cells
sampleCell <- function(obj, posCell, method = c("knn", "random"),
                       bgCellnumber=NULL){
  if(method=="knn"){
    if(is.null(obj$graphs$PCA)){
      stop("please run knn first")
    }
  }

  if(is.null(obj[["pmat"]])){
    stop("please add peak matrix to pagoda object by running addAtacObj")
  }

  ncell = nrow(obj[["pmat"]])
  ncell.pos = length(posCell)
  ncell.neg = min(length(posCell), ncell - ncell.pos)

  if(method == "knn"){
    # transform to adjcent matrix
    adj <- as_adjacency_matrix(obj$graphs$PCA)

    neg = setdiff(seq(ncell), posCell)
    dx = order(Matrix::colSums(adj[posCell, neg]),
               decreasing = TRUE)[1:ncell.neg]
    neg = neg[dx]
  }
  else if(method == "random") {
    neg = setdiff(seq(ncell), posCell)
    if(!is.null(bgCellnumber)){
      sampleSize = bgCellnumber
    } else{
      sampleSize = min(length(posCell), length(neg))
    }
    neg = sample(neg, sampleSize)
  }
  return(neg)
}

#' identify differential accessible region based on pagoda object
#' @param obj pagoda object
#' @param clusters factor contain cluster information for each cell
#' @param bgCellnumer background cell number
#' @export
#' @import edgeR
idenDAR <- function(obj, clusters, neg=NULL, method=c("knn", "random"), bcv=0.1,
  test.method=c("exactTest", "lrt", "qlf"), n.cores=2, bgCellnumber=NULL){
  require(edgeR)

  if(is.null(obj[["pmat"]])){
    stop("please add peak matrix to pagoda object")
  }

  if(is.null(clusters)){
    stop("please provide cluster vector")
  }

  pmat <- obj[["pmat"]]
  ncell <- nrow(pmat)

  ADRtest <- function(targetCL, pmat, method=c("knn", "random"), test.method=c("exactTest", "lrt", "qlf"), bcv=0.1){
    pos.index <- which(clusters == targetCL)
    if(method=="knn"){
      neg.index <- sampleCell(obj, posCell=pos.index, method=method, bgCellnumber=bgCellnumber)
    } else{
      neg.index <- sampleCell(obj, posCell=pos.index, method="random", bgCellnumber=bgCellnumber)
    }

    ncell.pos = length(pos.index)
    ncell.neg = length(neg.index)

    pmat.pos = pmat[pos.index, ,drop=FALSE]
    pmat.neg = pmat[neg.index, ,drop=FALSE]

    # performing test
    data = data.frame(Matrix::colSums(pmat.neg), Matrix::colSums(pmat.pos))
    group <- factor(c(1,2))
    design <- model.matrix(~group)
    y <- edgeR::DGEList(counts=data, group=group)

    if(test.method == "lrt"){
      fit <- glmFit(y, design, dispersion=bcv^2);
      diffADR <- glmLRT(fit,coef=2)$table;
    }else if(test.method == "qlf"){
      diffADR <- glmQLFit(y, design, dispersion=bcv^2)$table;
    }else{
      diffADR <- exactTest(y, dispersion=bcv^2)$table;
    }
    # sorting
    diffADR <- diffADR[order(diffADR$PValue), ]
    diffADR$adjPValue <- p.adjust(diffADR$PValue, method="BH")
    return(diffADR)
  }

  nclusters <- length(levels(clusters))
  # testing
  ADR <- parallel::mclapply(1:nclusters, function(r){
    targetCL <- r
    ADRtest(targetCL, pmat, method=method, test.method=test.method)
  }, mc.cores=n.cores)

  names(ADR) <- levels(clusters)
  obj[["DAR"]] <- ADR
  return(obj)
}



