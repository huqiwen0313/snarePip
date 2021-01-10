# basic functions for data processing and transformation

#' return sample library id based on experiment ID
getLibID <- function(linkTable, sampleName, assayfix=FALSE){
  
  if(assayfix){
    RNAtable <- linkTable[grep("SPL-R|SNARE2-R", linkTable$Experiment_ID), ]
    atacTable <- linkTable[grep("SPL-AC|SNARE2-AC", linkTable$Experiment_ID), ]
    # RNA table
    RNAtable$Experiment_ID <- gsub("[-|_]SPL-R", "", RNAtable$Experiment_ID)
    RNAtable$Experiment_ID <- gsub("[-|_]SNARE2-R", "", RNAtable$Experiment_ID)
    RNAtable$Experiment_ID <- gsub("_P\\d+", "", RNAtable$Experiment_ID)
    
    # ATAC table
    atacTable$Experiment_ID <- gsub("-SPL-AC", "", atacTable$Experiment_ID)
    atacTable$Experiment_ID <- gsub("_SNARE2-AC", "", atacTable$Experiment_ID)
    atacTable$Experiment_ID <- gsub("-SNARE2-AC", "", atacTable$Experiment_ID)
    atacTable$Experiment_ID <- gsub("_P[0-9]*", "", atacTable$Experiment_ID)
    
    linkSampleTable <- rbind(RNAtable, atacTable)
  } else{
    linkSampleTable <- linkTable
  }
  
  
  #
  sampleName <- gsub("-SPL.*Ad1", "_", sampleName)
  sampleName <- gsub("\\.P", "_P", sampleName)
  sampleName <- gsub(".SPLiT", "", sampleName)
  sampleName <- gsub("_S.*", "", sampleName)
  sampleName <- gsub("_S.*", "", sampleName)
  sampleName <- gsub(".Ad1", "", sampleName)
  sampleName <- gsub("_P[0-9]*", "", sampleName)
  libID <- linkSampleTable[linkSampleTable$Experiment_ID==sampleName, ]$Library_ID
  return(libID)
}

getTissue <- function(linkTable, sampleName=NULL, libID=NULL, assayfix=FALSE){
  
  if(assayfix){
    RNAtable <- linkTable[grep("SPL-R|SNARE2-R", linkTable$Experiment_ID), ]
    atacTable <- linkTable[grep("SPL-AC|SNARE2-AC", linkTable$Experiment_ID), ]
    # RNA table
    RNAtable$Experiment_ID <- gsub("[-|_]SPL-R", "", RNAtable$Experiment_ID)
    RNAtable$Experiment_ID <- gsub("[-|_]SNARE2-R", "", RNAtable$Experiment_ID)
    RNAtable$Experiment_ID <- gsub("_P\\d+", "", RNAtable$Experiment_ID)
    
    # ATAC table
    atacTable$Experiment_ID <- gsub("-SPL-AC", "", atacTable$Experiment_ID)
    atacTable$Experiment_ID <- gsub("_SNARE2-AC", "", atacTable$Experiment_ID)
    atacTable$Experiment_ID <- gsub("-SNARE2-AC", "", atacTable$Experiment_ID)
    atacTable$Experiment_ID <- gsub("_P[0-9]*", "", atacTable$Experiment_ID)
    
    linkSampleTable <- rbind(RNAtable, atacTable)
  } else{
    linkSampleTable <- linkTable
  }
  
  if(!is.null(sampleName)){
    tissue <- unique(linkSampleTable[linkSampleTable$Experiment_ID_Short==sampleName, ]$Tissue)
    return(tissue)
  }
  if(!is.null(libID)){
    tissue <- unique(linkSampleTable[linkSampleTable$Library_ID==libID, ]$Tissue)
    return(tissue)
  }
}

#' merge sparse matrix
merge.sparse <- function(listMatrixes) {
  # takes a list of sparse matrixes with different columns and adds them row wise
  allRownames <- sort(unique(unlist(lapply(listMatrixes,rownames))))
  allColnames <- as.character(unlist(lapply(listMatrixes,colnames)))
  for (currentMatrix in listMatrixes) {
    newRowLocations <- match(rownames(currentMatrix),allRownames)
    indexes <- Matrix::which(currentMatrix>0, arr.ind = T)
    newRows <- newRowLocations[indexes[,1]]
    columns <- indexes[,2]
    newMatrix <- sparseMatrix(i=newRows,j=columns, x=currentMatrix@x,
                              dims=c(length(allRownames),max(columns)))
    if (!exists("matrixToReturn")) {
      matrixToReturn <- newMatrix
    }
    else {
      matrixToReturn <- cbind2(matrixToReturn,newMatrix)
    }
  }
  rownames(matrixToReturn) <- allRownames
  colnames(matrixToReturn) <- allColnames
  matrixToReturn
}


p2proc <- function(cd, n.cores=20, min.cells.per.gene=200, min.cell.size=200, n.odgenes=2e3, 
                   get.largevis=FALSE, make.geneknn=FALSE, perplexity=50, log.scale = TRUE, 
                   batch = NULL, nPcs = 30, k = 40, get.tsne = TRUE){
    counts <- gene.vs.molecule.cell.filter(cd, min.cell.size=min.cell.size, plot=F)
    rownames(counts) <- make.unique(rownames(counts))
    p2 <- Pagoda2$new(counts, n.cores = n.cores, batch=batch,
                      log.scale = log.scale, min.cells.per.gene=min.cells.per.gene)
    p2$adjustVariance(plot=F, gam.k = 10)
    p2$calculatePcaReduction(nPcs = nPcs, n.odgenes = n.odgenes,
                             maxit = 1000)
    p2$makeKnnGraph(k = k, type = "PCA", center = TRUE, weight.type = "none",
                    n.cores = n.cores, distance = "cosine")
    #p2$getKnnClusters(method = igraph::infomap.community, type = "PCA",
    #                  name = "infomap")
    p2$getKnnClusters(method = igraph::multilevel.community,
                      type = "PCA", name = "multilevel")
    if (get.largevis) {
      M <- 30
      p2$getEmbedding(type = "PCA", embeddingType = "largeVis",
                      M = M, perplexity = perplexity, gamma = 1/M, alpha = 1)
    }
    if (get.tsne) {
      if (perplexity > nrow(p2$counts)/5) {
        perplexity <- floor((nrow(p2$counts) - 1)/3)
        cat("perplexity is too large, reducing to", perplexity,
            "\\n")
      }
      p2$getEmbedding(type = "PCA", embeddingType = "tSNE",
                      perplexity = perplexity)
    }
    if (make.geneknn){
      p2$makeGeneKnnGraph()
    }
    return(p2)
}


plotCompositionBarplots <- function(groups, sample.factor=NULL, 
				    show.entropy=TRUE,show.size=TRUE, show.composition=TRUE,legend.height=0.2) {
  ## param checking
 if (is.null(groups)) {
    stop('groups factor on the cells needs to be specified')
  }
  
  groups <- as.factor(groups)
  #if(is.null(sample.factor)) {
  #  sample.factor <- conos.obj$getDatasetPerCell(); # assignment to samples
  #}
  
  xt <- table(sample.factor[match(names(groups),names(sample.factor))],groups)
  xt <- xt[rowSums(xt)>0,]; xt <- xt[,colSums(xt)>0]
  
  df <- reshape2::melt(xt); colnames(df) <- c("sample","cluster","f");  df$f <- df$f/colSums(xt)[as.character(df$cluster)]
  clp <- ggplot2::ggplot(df, ggplot2::aes(x=factor(cluster, levels=levels(groups)),y=f,fill=sample)) +
    ggplot2::geom_bar(stat='identity') + ggplot2::xlab('cluster') + ggplot2::ylab('fraction of cells') + ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(expand=c(0, 0))
  
  if(!show.size && !show.entropy)
    return(clp);
  
  # extract legend
  leg <- cowplot::get_legend(clp + ggplot2::theme(legend.position="bottom"))
  pl <- list(clp + ggplot2::theme(legend.position="none"));
  
  if(show.entropy) {
    if (!requireNamespace("entropy", quietly=T))
      stop("You need to install 'entropy' package to use 'show.entropy=T'")
    
    n.samples <- nrow(xt);
    ne <- 1-apply(xt, 2, entropy::KL.empirical, y2=rowSums(xt), unit=c('log2')) / log2(n.samples) # relative entropy
    enp <- ggplot2::ggplot(data.frame(cluster=factor(colnames(xt),levels=levels(groups)),entropy=ne), ggplot2::aes(cluster, entropy)) +
      ggplot2::geom_bar(stat='identity',fill='grey65') + ggplot2::ylim(0,1) +
      ggplot2::geom_hline(yintercept=1, linetype="dashed", color = "grey30") + ggplot2::theme_bw()
    pl <- c(pl,list(enp))
  }
  
  if(show.size) {
    szp <- ggplot2::ggplot(data.frame(cluster=factor(colnames(xt),levels=levels(groups)), cells=colSums(xt)), ggplot2::aes(cluster,cells)) +
      ggplot2::geom_bar(stat='identity') + ggplot2::scale_y_continuous(trans='log10') + ggplot2::theme_bw() + ggplot2::ylab('number of cells')
    pl <- c(pl,list(szp))
  }
  
  pp <- cowplot::plot_grid(plotlist=pl,ncol=1,rel_heights=c(1,rep(0.3,length(pl)-1)))
  pp2 <- cowplot::plot_grid(leg,pp,ncol=1,rel_heights=c(legend.height,1))
  
  return(pp2)
}


clusterComposition <- function(groups, sample.factor=NULL, show.entropy=TRUE,show.size=TRUE, show.composition=TRUE,legend.height=0.2) {
  ## param checking
  if (is.null(groups)) {
    stop('groups factor on the cells needs to be specified')
  }
  
  groups <- as.factor(groups)

  xt <- table(sample.factor[match(names(groups),names(sample.factor))],groups)
  xt <- xt[rowSums(xt)>0,]; xt <- xt[,colSums(xt)>0]
  df <- reshape2::melt(xt); colnames(df) <- c("sample","cluster","f")  
  df$f <- df$f/colSums(xt)[as.character(df$cluster)]
  return(df)
}


GetScrubletScores <- function(mat, min.molecules.per.gene=10, pythonPath="./",
                              method=c("scrublet", "doubletDetection")) {
  tf.in <- tempfile()
  tf.out <- tempfile()
  
  as_matrix <- function(mat){
     temp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
     prow <- mat@i+1
     pcol <- findInterval(seq(mat@x)-1, mat@p[-1])+1
     val <- mat@x
     
     for (i in seq_along(val)){
      temp[prow[i], pcol[i]] <- val[i]
    }
    
    row.names(temp) <- mat@Dimnames[[1]]
    colnames(temp) <- mat@Dimnames[[2]]
    return(temp)
  }
  dt <- mat[Matrix::rowSums(mat)>=min.molecules.per.gene,] %>% Matrix::t() %>% as_matrix() %>% data.table::data.table()
  data.table::fwrite(dt, file=tf.in)
  
  if(method=="scrublet"){
    cmd <- paste0(pythonPath, " -c 'import sys; import pandas; import scrublet; ",
                  "df = pandas.read_csv(\"", tf.in, "\"); scrub = scrublet.Scrublet(df); ",
                  "doublet_scores, predicted_doublets = scrub.scrub_doublets();",
                  "pandas.DataFrame(dict(score=doublet_scores, is_doublet=predicted_doublets)).to_csv(\"",tf.out,"\");'",
                  sep='')
    system(cmd, intern=F)
    x <- data.table::fread(tf.out,sep=',')[,2:3] %>% as.data.frame() %>% as.list() %>%
      lapply(`names<-`, colnames(mat))
  } else if(method == "doubletDetection"){
    cmd <- paste0(pythonPath, " -c 'import sys; import pandas; import doubletdetection; ",
                  "df = pandas.read_csv(\"", tf.in, "\"); clf = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=False, standard_scaling=True); ",
                  "doublets = clf.fit(df).predict(p_thresh=1e-7, voter_thresh=0.8);",
                  #"doublet_scores, predicted_doublets = scrub.scrub_doublets();",
                  "pandas.DataFrame(doublets).to_csv(\"",tf.out,"\");'",
                  sep='')
    system(cmd, intern=F)
    x <- data.table::fread(tf.out,sep=',')[-1, 2]
    x <- data.frame(cells=colnames(mat), doublet=x)
  }
  
  return(x)
}


# from https://github.com/khodosevichlab/Epilepsy19/blob/master/R/wrappers.R
GetPagodaWebApp <- function(p2, clusters, organism=NULL, additional.metadata=list(), 
                            verbose=T, go.sets=NULL, go.env=NULL, test.pathways=TRUE) {
  
  ExtractGoSets <- function(go.env) {
    names(go.env) %>% setNames(., .) %>% sccore:::plapply(function(x)
      list(properties = list(locked=T, genesetname=x, shortdescription=GO.db::GOTERM[[x]]@Term),
           genes = c(go.env[[x]])))
  }
  
  if (is.null(go.env)) {
    if (is.null(organism))
      stop("Either organism or go.env must be provided")
    
    if (verbose) cat("Generate go environment\n")
    go.env <- pagoda2::p2.generate.go(p2, organism=organism)
  }
  
  if (is.null(go.sets)) {
    if (verbose) cat("Generate genesets\n")
    
    
    go.sets <- ExtractGoSets(go.env)
  }
  
  if (verbose) cat("Generate de geneset\n")
  de.sets <- pagoda2::get.de.geneset(p2, groups = clusters, prefix = 'de_')
  
  if (test.pathways) {
    if (verbose) cat("Test pathway overdispersion\n")
    p2$testPathwayOverdispersion(setenv = go.env, verbose = verbose, correlation.distance.threshold = 0.8,
                                 recalculate.pca = F, min.pathway.size = 50, max.pathway.size = 1000)
  }
  
  if (verbose) cat("Create app\n")
  p2.web <- p2 %>%
    pagoda2::make.p2.app(
      dendrogramCellGroups = as.factor(clusters),
      geneSets = c(go.sets, de.sets),
      additionalMetadata=additional.metadata,
      show.clusters = T);
  
  if (verbose) cat("All done!\n")
  
  return(p2.web)
}

GetGenesetFraction_overall <- function (count.matrix, genes){
  umi.counts <- sort(Matrix::colSums(count.matrix), decreasing = T)
  presented.mit.genes <- intersect(genes, rownames(count.matrix))
  genes.frac <- sum(Matrix::colSums(count.matrix[presented.mit.genes, 
                                                 names(umi.counts)]))/sum(umi.counts)
  return(genes.frac)
}

