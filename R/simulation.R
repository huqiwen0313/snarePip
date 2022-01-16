# simulation related functions

#' calculate softmax
softmax <- function(v){
  sm <- exp(v)/sum(exp(v))
  return(sm)
}


#' calculate change for the 1st celltype
#' @param counts vector of counts for Refs
#' @param clevel abundance change between groups
#' @param n.cells number of cells
muPrime <- function(counts, clevel, n.cells){
  n.celltype <- length(counts)
  b <- log(counts/n.cells)
  
  counts[1] <- counts[1] + clevel
  counts[2:n.celltype] <- round((n.cells-counts[1])/(n.cells-1), 3)
  
  b_prime <- log(counts/n.cells)
  
  w <- rep(0, n.celltype)
  w[1] <- round((b_prime - b)[1], 3)
  return(list(w=w, b=b))
}


#' generate effect matrix
#' @param n_d number of covariant affect each celltype
#' @param n_k number of celltype affect by each covariant
generateEffectMatrix <- function(d=1^2, n.cell.types, n_d, n_k, rc=c(0.3, 0.5, 1)){
  de <- sample(1:d, n_d, replace=F)
  ke <- sample(1:n.cell.types, n_k, replace=F)
  
  w <- lapply(seq(1, d*n.cell.types, 1), function(r){
    sample(rc, 1)
  }) %>% unlist() %>% matrix(., nrow=d, ncol=n.cell.types)
  
  return(w)
}


#' generate simulated composition matrix based on multinomial distribution (scCODA)
#' @param n.cell.types number of celltypes in each group
#' @param ncells number of cells for each sample
#' @param nsamples vector of number of samples c(5, 5)
#' @param b bias coefficients
#' @param w Effect matrix
generateSimMatrix <- function(n.cell.types, ncells, nsamples, ngroups=1, sigma=NULL,
                            b=NULL, w=NULL){
  if(is.null(b)){
    b = runif(n.cell.types, -3, 3)
  }
  
  if(is.null(w)){
    n_d <- sample(1:ngroups-1, 1)
    n_k <- sample(1:n.cell.types-1, 1)
    w <- generateEffectMatrix(d=ngroups^2, n.cell.types, n_d, n_k)
  }
  
  if(is.null(sigma)){
    sigma <- diag(n.cell.types)*0.05
  }
  
  x<- matrix(0L, nrow=sum(nsamples), ncol=ngroups)
  x[(nsamples[1]+1):sum(nsamples), ] <- 1
    
  #y <- matrix(0L, nrow=sum(nsamples), ncol=n.cell.types)
  
  y <- lapply(1:(ngroups+1), function(r){
    yc <- lapply(1:nsamples[r], function(j){
      if(r==1){
        c <- 0
      } else{
        c <- nsamples[r]
      }
      alpha <- mvtnorm::rmvnorm(1, mean=t(x[c+j, ]) %*% w+b, sigma=sigma)
      z = rmultinom(n=1, size=ncells, prob=softmax(alpha))
      t(z)
    }) %>% unlist() %>% matrix(., ncol=n.cell.types, byrow=T)
    #c <- c + nsamples[r]
    yc
  }) %>% do.call(rbind, .)
  
  return(list(y=y, x=x, w=w, b=b))
}

#' generate simulated composition datasets
#' @param base mean abundance for the 1st celltype
#' @param clevel abundance change between groups
#' @param re.ratio fraction of celltypes that have region effect
#' @param sample.re.ratio ratio of sample that have region effect
#' @param bias fraction of samples that have region effect between disease and 
#'            control - 0.5 means equal fraction
#' @param c.selected number vector showing which celltype to select
#' @param n.rep number of replicates for each parameter
#' @param batch.sime TRUE/FALSE, if simulate continuous batch effect for specific celltype
#' @param batch.factor vector to specify batch factor for each sample
#' @param batch.model pass model to predict sample cell fraction (X) based on batch factor (Y)
#' @return list of simulated count matrix with different replications
#'          and correspondent paramters  
generateSimData <- function(n.celltypes, n.cells, n.samples, base, clevel,
                            re.ratio=NULL, sample.re.ratio=0.2, bias=0.5,
                            c.selected=NULL, n.rep=1, batch.sim=FALSE,
                            batch.factor=NULL, batch.model=NULL){
  # generate counts
  b <- rep((n.cells - base)/(n.celltypes-1), n.celltypes)
  b[1] <- base
  b.t <- round(log(b/n.cells), 3)
  w <- muPrime(b, clevel, n.cells)$w
  
  lf = round(log2((clevel+base)/base), 2)
  
  # generate simulated matrix
  sigma <- diag(n.celltypes) * 0.05
  
  # region effect
  if(!is.null(re.ratio)){
    # number of celltypes that have region effect
    nc <- ceiling(sum(n.celltypes) * re.ratio)
    # number of samples that have region effect
    ns <- ceiling(sum(n.samples) * sample.re.ratio)
    
    # random sample celltype if NULL
    if(is.null(c.selected)){
      celltype <- seq(1, n.celltypes, 1)[sample(1:n.celltypes, nc, replace=F)]
    } else{
      celltype <- seq(1, n.celltypes, 1)[c.selected]
    }
    
    # get sample information
    n.control <- ceiling(ns*(1-bias))
    n.case <- ceiling(ns*bias)
    
    if(n.control > 0){
      if(n.control >= n.samples[1]){
        control.index <- 1:n.samples[1]
      } else{
        control.index <- sample(1:n.samples[1], size=n.control, replace=F)
      }
    } else{
      control.index <- 0
    }
    
    if(n.case >0){
      if(n.case >= n.samples[1]){
        case.index <- (n.samples[1]+1):sum(n.samples)
      } else{
        case.index <- sample((n.samples[1]+1):sum(n.samples), size=n.case, replace=F)
      }
    } else{
      case.index <- 0
    }
    sample.index <- unique(c(control.index, case.index)) %>% .[.>0]
    
  }
  
  mats <- lapply(1:n.rep, function(r){
    mat <- generateSimMatrix(n.cell.types=n.celltypes, ncells=n.cells, nsamples=n.samples,
                             sigma=sigma, b=b.t, w=w)
    
    rownames(mat$y) <- paste("p", seq(1, nrow(mat$y), 1), sep="")
    colnames(mat$y) <- paste("celltype", seq(1, ncol(mat$y)), sep="")
    
    if(!is.null(re.ratio)){
      mat$y[sample.index, celltype] <- 0
    }
    
    if(batch.sim){
      cell.frac <- predict(model, data.frame(X=batch.factor))
      # norm
      cell.frac <- (cell.frac - min(cell.frac))/sum(cell.frac - min(cell.frac))
      trt.cell.frac <- cell.frac[1:n.samples[1]]/sum(cell.frac[1:n.samples[1]])
      control.cell.frac <- cell.frac[(n.samples[1]+1):sum(n.samples)]/sum(cell.frac[(n.samples[1]+1):sum(n.samples)])
      
      if(is.null(c.selected)){
        c.selected <- seq(1, n.celltypes, 1)[sample(1:n.celltypes, nc, replace=F)]
      }
      
      # re-assign cell count
      for(i in c.selected){
        total.trt <- sum(mat$y[1:n.samples[1], i])
        mat$y[1:n.samples[1], i] <- round(total.trt*trt.cell.frac, 0)
          
        total.control <- sum(mat$y[(n.samples[1]+1):sum(n.samples), i])
        mat$y[(n.samples[1]+1):sum(n.samples), i] <- round(total.control*control.cell.frac, 0)
      }
       
    }
    
    mat$y
  })
  names(mats) <- paste("rep", 1:n.rep, sep="_")
  
  d.groups <- c(rep(FALSE, n.samples[1]), rep(TRUE, n.samples[2]))
  names(d.groups) <- rownames(mats[[1]])
  
  sample.groups <- as.factor(setNames(c(rep("Ref", n.samples[1]), rep("Disease", n.samples[2])), rownames(mats[[1]])))
  if(!is.null(re.ratio)){
    params <- list(n.celltypes=n.celltypes, n.cells=n.cells, n.samples=n.samples, 
                  base=base, clevel=clevel, lf=lf, b=b.t, w=w, d.groups=d.groups,
                  sample.groups=sample.groups, sample.index=sample.index,
                  rcelltypes.index=celltype, re.ratio=re.ratio, sample.re.ratio=sample.re.ratio,
                  bias=bias)
  } else{
    params <- list(n.celltypes=n.celltypes, n.cells=n.cells, n.samples=n.samples, 
                  base=base, clevel=clevel, lf=lf, b=b.t, w=w, d.groups=d.groups,
                  sample.groups=sample.groups, re.ratio=re.ratio, sample.re.ratio=sample.re.ratio,
                  bias=bias)
  }
  
  return(list(count=mats, params=params))
}

#' plot celltype proportion for simulated data
#' @param mat count matrix
plotSimProportions <- function(mat, sample.groups, legend.position = "right", notch=FALSE, alpha=0.1, palette=NULL,
                               show.significance=TRUE) {
  
  df <- as.data.frame(mat/rowSums(mat)) 
  #df$group <- mat$params$sample.groups
  df$group <- sample.groups
  df.melt <- reshape2::melt(df, id.vars="group")
  
  gg <- ggplot(df.melt, aes(x=variable, y=value, by=group)) +
    geom_boxplot(position=position_dodge(), outlier.shape = NA, notch=notch) +
    ylab("% cells per sample") +
    xlab("") +
    theme_bw() +
    #theme_legend_position(legend.position) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
          legend.title=element_blank()) +
    geom_point(position=position_jitterdodge(jitter.width=0.15), aes(col=group), alpha=alpha)
    #scale_y_continuous(expand=c(0, 0), limits=c(0, (max(df.melt$value) + 1)))
  
  if(show.significance) gg <- gg + ggpubr::stat_compare_means(aes(group = group), label = "p.signif")  # willcox test
  if(!is.null(palette)) gg <- gg+ scale_color_manual(values=palette)
  gg
}

#' calculate performance metrics
#' @param df data.frame with each column records tp, tn, pf, pn
#' @return df with metric added
getMetric <- function(df){
  df <- lapply(1:nrow(df), function(r){
    data <- df[r, ]
    tpr = (data$tp / (data$tp + data$fn))
    tnr = (data$tn / (data$tn + data$fp))
    precision = (data$tp / (data$tp + data$fp))
    fdr = (data$fp / (data$tp + data$fp))
    acc = (data$tp + data$tn) / (data$tp + data$tn + data$fp + data$fn)
    youdern = tpr + tnr - 1
    f1.score = 2 * (tpr * precision / (tpr + precision))
    mcc <- ((data$tp * data$tn) - (data$fp * data$fn)) / sqrt((data$tp + data$fp) * (data$tp + data$fn) * (data$tn + data$fp) * (data$tn + data$fn))
    
    data$mcc <- mcc
    data$tpr <- tpr
    data$tnr <- tnr
    data$precision <- precision
    data$fdr <- fdr
    data$acc <- acc
    data$youdern <- youdern
    data$f1.score <- f1.score
    
    data
  }) %>% dplyr::bind_rows() %>% as.matrix()
  df[which(is.na(df))] <- 0
  return(as.data.frame(df))
}


#' caculate performance for simulated dataset
#' @param mat simulation object generated by generateSimData function
simPerformance <- function(mat, n.cell.counts=1000, n.iter=1000, p.cutoff=0.001, ref.celltype=NULL){
  n.rep <- length(mat$count)
  
  df <- lapply(1:n.rep, function(r){
    t <- resampleContrast(mat$count[[r]], mat$params$d.groups, n.cell.counts=1000, n.iter=1000, ref.celltype=ref.celltype)
    pvals <- getCellSignificance(t$balances, ref.celltype=ref.celltype)
    sigCelltypes <- names(pvals)[pvals<p.cutoff]
    
    ide <- factor(names(pvals) %in% sigCelltypes, levels=c("TRUE", "FALSE"))
    true <- factor(c(TRUE, rep(FALSE, (mat$params$n.celltypes-1))), levels=c("TRUE", "FALSE"))
    confusion <- table(ide, true)
    
    tp <- confusion[1]
    fp <- confusion[3]
    fn <- confusion[2]
    tn <- confusion[4]
    #mcc <- mltools::mcc(true, ide)
    data.frame(rep=r, tp=tp, fp=fp, fn=fn, tn=tn)
  }) %>% dplyr::bind_rows()
  
  df$ncell <- mat$params$n.cells
  df$ncelltype <- mat$params$n.celltypes
  df$nsamples <- sum(mat$params$n.samples)
  df$mu <- mat$params$base
  df$mu.pr <- mat$params$clevel
  df$lf <- mat$params$lf
  
  df <- getMetric(df)
  return(df)
}
