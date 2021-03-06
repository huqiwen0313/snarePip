% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motifs.R
\name{enrichMotifsDAR}
\alias{enrichMotifsDAR}
\title{identify overrepresentative motifs in DAR regions for each cluster}
\usage{
enrichMotifsDAR(
  obj,
  sampleMethod = "match.GC",
  PvalueCutoff = 0.05,
  peakCutoff = 500,
  nsample = 10000,
  genome = "BSgenome.Hsapiens.NCBI.GRCh38",
  use.pvalues = FALSE
)
}
\arguments{
\item{obj}{pagoda2 object}

\item{PvalueCutoff}{pvalue cutoff for significant DAR regions}
}
\description{
identify overrepresentative motifs in DAR regions for each cluster
}
