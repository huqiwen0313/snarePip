#' simulate doublet matrix
simulateDoubletMatrix <- function(mat, sampleRatio1=0.5, sampleRatio2=0.5, nSample = 1000,
                                  clusters=NULL, sampleCluster=TRUE, clusterRatio=0.25){
  sampleMatrix <- function(mat, sampleRatio = 0.5){
    total <- length(mat@x)
    sampledloc <- floor(total * (1-sampleRatio))
    mat@x[sample(seq_len(total), sampledloc)] <- 0
    mat <- drop0(mat)
    mat
  }
  if(is.null(clusters) & sampleCluster==TRUE){
    stop("please provide cluster information")
  }

  if(sampleCluster){
    clusterName <- unique(clusters)

    simulatedMat <- lapply(1:length(clusterName), function(r){
      cluster1Cell <- names(clusters[clusters==r])
      pairedCluster <- sample(clusterName[-which(clusterName == r)], 1)
      cluster2Cell <- names(clusters[clusters==pairedCluster])
      nSample <- max(floor(length(cluster1Cell)*clusterRatio), floor(length(cluster2Cell)*clusterRatio))
      mat1 <- mat[rownames(mat) %in% cluster1Cell, ]
      mat2 <- mat[rownames(mat) %in% cluster2Cell, ]

      id1 <- sample(seq_len(nrow(mat1)), nSample, replace = TRUE)
      id2 <- sample(seq_len(nrow(mat2)), nSample, replace = TRUE)
      simulatedMat <- sampleMatrix(mat1[id1, ], sampleRatio = sampleRatio1) +
        sampleMatrix(mat2[id2, ], sampleRatio = sampleRatio2)
      rownames(simulatedMat) <-  make.unique(paste("doublet", rownames(simulatedMat), sep="."))
      simulatedMat
    }) %>% Reduce("rbind", .)

  } else{
    id1 <- sample(seq_len(nrow(mat)), nSample, replace = TRUE)
    id2 <- sample(seq_len(nrow(mat)), nSample, replace = TRUE)
    simulatedMat <- (sampleMatrix(mat[id1, ], sampleRatio = sampleRatio1) +
      sampleMatrix(mat[id2, ], sampleRatio = sampleRatio2))/2
    rownames(simulatedMat) <-  make.unique(paste("doublet", rownames(simulatedMat), sep="."))
  }
  return(simulatedMat)
}

computeKnn <- function(data, query=NULL, k=50, includeSelf=FALSE){
  require("nabor")

  if(is.null(query)){
    query <- data
    Self <- TRUE
  } else{
    Self <- FALSE
  }

  if(Self & !includeSelf){
    knnIdx <- nabor::knn(data, query, k=k + 1)$nn.idx
    knnIdx <- knnIdx[, -1, drop = FALSE]
  } else{
    knnIdx <- nabor::knn(data, query, k=k)$nn.idx
  }
  return(knnIdx)
}

#' calculate score for potential doublet
#' @export
calDoubletScore <- function(obj, plot=T, plotDoubleEmbedding=T, plotDepthDistr=T, k=30){
  require(pagoda2)
  if(is.null(obj[["pmat"]])){
    stop("please provide peak matrix")
  }
  pmat <- obj[["pmat"]]

  # get original clusters
  obj <- RunJD(obj)
  obj <- normJD(obj)
  obj$counts <- Matrix::Matrix(obj[["jmat"]]$jmat, sparse=TRUE)
  obj$adjustVariance(plot=F, gam.k=10)
  obj$counts <- Matrix::Matrix(obj[["njmat"]], sparse=TRUE)
  obj$misc$odgenes <- colnames(obj[["njmat"]])
  obj$calculatePcaReduction(nPcs=50, n.odgenes=3e6)
  obj$makeKnnGraph(k=15, type='PCA', center=T, distance='cosine')
  p2graph <- obj$graphs$PCA
  clusters <- igraph::membership(conos::leiden.community(p2graph, resolution=2))
  obj$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=20,verbose=F)

  # simulate doublet
  nSample <- floor(nrow(pmat) * 0.25)
  simulatedMat <- simulateDoubletMatrix(pmat, sampleRatio1=0.5, sampleRatio2=0.5,
                                        nSample=nSample, sampleCluster=TRUE, clusters=clusters)

  p2Mixed <- Pagoda2$new(t(simulatedMat), log.scale=TRUE)
  p2Mixed[["pmat"]] <- rbind(simulatedMat)
  p2Mixed <- RunJD(p2Mixed)
  p2Mixed <- normJD(p2Mixed)
  p2Mixed$counts <- Matrix::Matrix(p2Mixed[["jmat"]]$jmat, sparse=TRUE)
  p2Mixed$adjustVariance(plot=F, gam.k=10)
  p2Mixed$counts <- Matrix::Matrix(p2Mixed[["njmat"]], sparse=TRUE)
  p2Mixed$misc$odgenes <- colnames(p2Mixed[["njmat"]])
  p2Mixed$calculatePcaReduction(nPcs=50, n.odgenes=3e6)
  p2Mixed$makeKnnGraph(k=15, type='PCA', center=T, distance='cosine')
  p2Mixed$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=20,verbose=F)

  # get embedding
  doubletEmbedding <- p2Mixed$embeddings$PCA$tSNE
  sigletEmbedding <- obj$embeddings$PCA$tSNE

  if(plot){
    doublet.cells <- rownames(simulatedMat)
    single.cells <- rownames(pmat)
    new_clusters <- setNames(c(rep("doublet", length(doublet.cells)), rep("single Cells", length(single.cells))),
                             c(doublet.cells, single.cells))

    p2Mixed$embeddings$PCA$tSNE <- rbind(p2Mixed$embeddings$PCA$tSNE, obj$embeddings$PCA$tSNE)

    par(mfrow=c(1, 2))
    obj$plotEmbedding(type='PCA',embeddingType='tSNE',groups=clusters,
                          show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,
                          mark.cluster.cex=1,alpha=0.5,main='original embedding, tSNE', quiet=T)
    p2Mixed$plotEmbedding(type='PCA',embeddingType='tSNE',groups=new_clusters,
                          show.legend=T,mark.clusters=F,min.group.size=1,shuffle.colors=F,
                          mark.cluster.cex=1,alpha=0.5,main='embedding with simulated doublet, tSNE', quiet=T)
  }

  knnPoints <- computeKnn(sigletEmbedding, doubletEmbedding, k)
  knnDistr <- table(as.vector(knnPoints))
  knnCount <- rep(0, nrow(sigletEmbedding))
  knnCount[as.integer(names(knnDistr))] <-  knnCount[as.integer(names(knnDistr))] + knnDistr

  scale <- 10000/length(knnCount)
  pvalBinom <- unlist(lapply(1:length(knnCount), function(r){
    countKnn <- round(knnCount[r] * scale)
    sumKnn <- round(sum(countKnn) * scale)
    pbinom(countKnn - 1, sumKnn, 1/10000, lower.tail = FALSE)
  }))

  padj <- p.adjust(pvalBinom, method = "bonferroni")
  doubletScore <- -log10(pmax(padj+0.0000001, 1e-500))
  doubletEnrich <- 10000*((knnCount/sum(knnCount)) / (1/nrow(sigletEmbedding)))/length(knnCount)

  douletenrichment <- setNames(doubletEnrich, rownames(sigletEmbedding))

  if(plotDoubleEmbedding){
    par(mfrow=c(1, 1))
    obj$plotEmbedding(type='PCA',embeddingType='tSNE',
                      colors=douletenrichment, shuffle.colors=F,mark.cluster.cex=1,alpha=0.5,
                      main="Doublet enrichment score")
  }

  #### temporary
  if(plotDepthDistr){
    # binarize
    bmat <- pmat[pmat>0] <- 1

    readDepth <- Matrix::rowMeans(pmat)
    score.depth <- data.frame(cell=names(douletenrichment), enrichmentScore=douletenrichment)
    score.depth <- merge(score.depth, data.frame(cell=names(readDepth), depth=readDepth), by=c(1))
    cluster.table <- data.frame(cell=names(clusters), cluster=as.numeric(clusters))
    score.depth <- merge(score.depth, cluster.table, by=c(1))

    par(mfrow=c(ceiling(nrow(table(score.depth$cluster))/3), 3))
    lapply(1:nrow(table(score.depth$cluster)), function(r){
      score.depth.cluster <- score.depth[score.depth$cluster==r, ]
      plot(score.depth.cluster$enrichmentScore, score.depth.cluster$depth, pch=16,
           xlab="enrichment score", ylab="average depth", main=paste("cluster", r, sep=" "))
    })
  }

  return(douletenrichment)
}
