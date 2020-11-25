#' Creat additional slots for atac-seq analysis extended from an existing object
#' @param obj existing object (paagoda2)
#' @param pmat peak matrix - rows are cells colums are peaks (chr:start-end)
#' @param frament fragment file that has the same format as
#'                 https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments
#' @param promoterRegion GRanges object containing the promoter regions
#' @param TSS GRanges object containing tss locations
#' @export
#' @import tidyr
addAtacObj <- function(obj, pmat=NULL, fragments=NULL, promoterRegion=NULL, TSS=NULL,
                       blacklist=NULL, exons=NULL, transcripts=NULL, enhancers=NULL, tool="pagoda2"){
  if(is.null(pmat) & is.null(fragments)){
    print("cell by peak matrix is or fragment file is NULL")
  }
  if(!is.null(pmat)){
    obj[["pmat"]] <- pmat
    peaks <- data.frame(chr = colnames(pmat)) %>%
      tidyr::separate(col = chr, into = c("chr", "start", "end"), sep = ":|-")
    obj[["peaks"]] <- GRanges(peaks[, 1], IRanges(as.numeric(peaks[, 2]), as.numeric(peaks[, 3])))
  } else{
    obj[["pmat"]] <- NULL
    obj[["peaks"]] <- NULL
  }

  obj[["fragments"]] <- fragments
  obj[["promoterRegion"]] <- promoterRegion
  obj[["TSS"]] <- TSS
  obj[["blacklist"]] <- blacklist
  obj[["exons"]] <- exons
  obj[["transcripts"]] <- transcripts
  obj[["enhancers"]] <- enhancers
  return(obj)
}

#' generate motif to feature(peak) matrix, row motifs, column matching score for a feature(peak)
#' @param obj pagoda2 object
#' @param outype types to count motifs (scores, matches)
#' @return obj with motif block
#' @export
addMotifMatrix <- function(obj, genome="BSgenome.Hsapiens.NCBI.GRCh38", outype="matches"){
  require(JASPAR2018)
  require(TFBSTools)
  require(motifmatchr)

  if(is.null(obj[["peaks"]])){
    stop("please import peak matrix to object using addAtacObj function")
  }
  peaks <- obj[["peaks"]]
  print("extrac peaks done...")
  if(genome == "BSgenome.Hsapiens.NCBI.GRCh38"){
    require(BSgenome.Hsapiens.NCBI.GRCh38)
    genomeSeq <- BSgenome.Hsapiens.NCBI.GRCh38
    seqnames(genomeSeq) <- paste("chr", GenomicRanges::seqnames(genomeSeq), sep="")
    #species code
    spcode <- 9606

    # remove uncharaterized chromosomes
    chrlist <- paste("chr", seq(1:22), sep="")
    chrlist <- c(chrlist, "chrMT", "chrX", "chrY")
    peaks <- peaks[as.character(peaks@seqnames) %in% chrlist]
  }

  pwms <- TFBSTools::getMatrixSet(JASPAR2018, opts = list(species = spcode, all_versions = FALSE))
  print("matching motifs to motif profiles....")
  motifMatch <- motifmatchr::matchMotifs(pwms = pwms, peaks,
                                       genome = genomeSeq, out = outype)
  if(outype=="scores"){
    motifMatrix <- motifmatchr::motifScores(motifMatch)
  } else if(outype == "matches"){
    motifMatrix <- motifmatchr::motifMatches(motifMatch)
  }
  rownames(motifMatrix) <- paste0(as.character(seqnames(peaks)), ":",
                                   start(peaks), "-", end(peaks))
  motifname <- unlist(lapply(pwms, function(r){
    r@name
  }))

  obj[["motifMatrix"]] <- list(matrix=motifMatrix, motifname=motifname)
  return(obj)
}
