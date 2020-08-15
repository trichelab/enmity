#' longest common (gapped) substring for a sequence of strings
#' 
#' This is built for comfort, not for speed.
#' Don't use it to build phylogenetic trees.
#'
#' @param ...   the strings
#' @param sep   their internal separator (default is "") 
#' 
#' @return      their LCS
#' 
#' @export
LCS <- function(..., sep="") {

  paste(Reduce(intersect, strsplit(unlist(...), sep)), collapse=sep)

}
