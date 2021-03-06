% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadBpFiles.R
\name{loadBpFiles}
\alias{loadBpFiles}
\title{load results of BiocParallel jobs when !is.na(bpresultdir(BPPARAM))}
\usage{
loadBpFiles(BPPARAM, simplify = TRUE)
}
\arguments{
\item{BPPARAM}{a subclass of SnowParam (e.g. SerialParam, MulticoreParam)}

\item{simplify}{if a list of equal sized arrays is found, cbind it? (TRUE)}
}
\value{
\preformatted{       whatever is in the result directory, tidied up if possible
}
}
\description{
load results of BiocParallel jobs when !is.na(bpresultdir(BPPARAM))
}
