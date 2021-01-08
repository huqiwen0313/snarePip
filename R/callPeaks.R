# peaks calling functions

#' call peaks based on macs2
#' obj: p2 object record the fragment file paths
#' @param macs.options: macs option for peak calling - adapted from snapATAC option
#'	(https://github.com/r3fang/SnapATAC/blob/master)
#' @export

callMACs <- function(fragmentFileDir, output.prefix="", path.to.macs, gsize,
                     buffer.size=500, out.dir=getwd(),
                     macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
                     keep.simple=TRUE){
  frag.dir <- fragmentFileDir
  #frag.dir <- obj[["fragmentsPaths"]]
  #if(is.null(frag.dir)){
  #  stop("error: fragment file doest exist - please push fragment file directory into p2 object")
  #}

  if(missing(path.to.macs)){
    stop("please provide macs runing path")
  }

  if(missing(gsize)){
    stop("gsize is missing");
  }

  # run macs2
  message("start call peaks using macs2......")
  flag <- system2(command=path.to.macs,
                  args=c("callpeak",
                        "-t", frag.dir,
                        "-f", "BED",
                        "-g", gsize,
                        macs.options,
                        "-n", output.prefix,
                        "--outdir", out.dir))
  if (flag != 0) {
    stop("'MACS' call failed");
  }

  if(keep.simple){
    system(paste("rm ", out.dir, "/", output.prefix, "_control_lambda.bdg", sep=""))
    system(paste("rm ", out.dir, "/", output.prefix, "_peaks.xls", sep=""));
    system(paste("rm ", out.dir, "/", output.prefix, "_summits.bed", sep=""));
  }

  peak.file <- paste(out.dir, "/", output.prefix, "_peaks.narrowPeak", sep="")
  if(!file.size(peak.file) == 0){
	  peaks <- read.table(paste(out.dir, "/", output.prefix, "_peaks.narrowPeak", sep=""))
	  return(peaks)
  } else{
	  return(data.frame())
  }
}

#' call peaks based on assigned cell clusters
#' @param obj p2 object
#' @param clusters cell annotation
#' @export
callMACsCluster <- function(obj, clusters=NULL, path.to.macs, gsize,
                            buffer.size=500, out.dir=getwd(),
                            macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
                            ncores=1, keep.simple=TRUE,
                            min.cells=10, includeOthers=TRUE){
  require(data.table)
  if(is.null(clusters)){
    stop("please provide cell clusters")
  }

  if(missing(path.to.macs)){
    stop("please provide macs runing path")
  }

  if(missing(gsize)){
    stop("gsize is missing");
  }
  frag.file <- obj[["fragmentsPaths"]]
  frag.dir <- strsplit(frag.file, split="/")[[1]]
  sampleName <- gsub(".fragments.bed", "", frag.dir[length(frag.dir)])
  frag.dir <- paste(frag.dir[-c(length(frag.dir))], collapse='/')

  if(is.null(frag.dir)){
    stop("error: fragment file doest exist - please push fragment file directory into p2 object")
  }
  # read fragment file
  #fragments <- read.table(frag.dir, sep="\t")
  fragments <- data.table::fread(frag.file, showProgress = FALSE)

  # extract cells
  allcells <- unique(fragments$V4)
  if(includeOthers){
    otherCells <- allcells[-1*which(allcells %in% names(clusters))]
    orderCells <- c(names(clusters), otherCells)
    clusters <- c(as.numeric(clusters), rep(max(as.numeric(clusters))+1, length(otherCells)))
    clusters <- setNames(clusters, orderCells)
    #clusters <- as.factor(clusters)
  }

  # call peaks
  clusterName <- names(table(clusters))
  peaks.gr.all <-  parallel::mclapply(1:length(clusterName), function(r){
    print(paste("calling peaks in cluster ", r, " ......", sep=""))
    clusterCells  <- names(clusters[clusters==clusterName[r]])
    num.cells = length(clusterCells)

    if(num.cells < min.cells){
      return(GenomicRanges::GRanges())
    }

    # extract fragments in cluster
    clusterBed <- fragments[fragments$V4 %in% clusterCells, ]
    # write to directory
    cluster.frag.dir <- file.path(frag.dir, paste(sampleName, r, "fragments.bed", sep="."))
    data.table::fwrite(clusterBed,
                file=cluster.frag.dir,
                sep = "\t", row.names=FALSE, col.names=FALSE)

    cluster.peaks.df <- callMACs(
      cluster.frag.dir,
      output.prefix=paste(paste(sampleName, r, sep=".")),
      path.to.macs=path.to.macs,
      gsize=gsize,
      macs.options=macs.options,
      out.dir=out.dir,
      keep.simple=keep.simple
    )

    if(nrow(cluster.peaks.df) > 0L){
      peaks.gr <- GenomicRanges::GRanges(cluster.peaks.df[,1], IRanges(cluster.peaks.df[,2], cluster.peaks.df[,3]));
    }else{
      peaks.gr <- GenomicRanges::GRanges();
    }
    peaks.gr
  }, mc.cores=ncores)
  peaks.gr.all <- GenomicRanges::reduce(do.call(base::c, peaks.gr.all))

  # write merged peak files
  peaks.merged <- as.data.frame(peaks.gr.all)[, 1:3]
  data.table::fwrite(peaks.merged,
                     file=paste(file.path(out.dir, sampleName), "combined.peaks.bed", sep="."),
                     sep = "\t", row.names=FALSE, col.names=FALSE)
  return(peaks.gr.all)
}




