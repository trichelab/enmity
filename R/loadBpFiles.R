#' load results of BiocParallel jobs when !is.na(bpresultdir(BPPARAM))
#' 
#' @param   BPPARAM   a subclass of SnowParam (e.g. SerialParam, MulticoreParam)
#' @param   simplify  if a list of vectors is found, cbind them? (TRUE) 
#' 
#' @return            whatever is in the result directory, tidied up if possible
#' 
#' @import  BiocParallel
#' 
#' @export
loadBpFiles <- function(BPPARAM, simplify=TRUE) { 

  resultDir <- bpresultdir(BPPARAM)
  stopifnot(.checkResultDir(resultDir))
  files <- list.files(resultDir, patt="^BP.*Rda$")
  res <- lapply(files, .loadBpFile, path=resultDir, simplify=simplify)
  if (simplify) res <- .simplifyres(res)
  return(res) 

}


# load one BiocParallel saved result
.loadBpFile <- function(file, path=".", simplify=TRUE) {
  
  stopifnot(.checkResultDir(path))
  res <- get(load(file.path(path, file)))
  if (simplify) res <- .simplifyres(res)
  return(res) 

}


# common case 
.simplifyres <- function(res) {
 
  if (.equalrows(res) & .equalcols(res)) {
    if (.hasnames(res)) nameses <- .getnames(res)
    res <- do.call(cbind, res)
    if (exists("nameses")) colnames(res) <- nameses
  }
  
  return(res)

}


# see if a list value has name attributes for cbinding etc. 
.attrnames <- function(x) names(attributes(x))


# see if an object has an attribute 
.hasattr <- function(x, y="name") return(y %in% .attrnames(x))


# common case
.hasnames <- function(x) all(sapply(x, .hasattr, "name"))


# common case
.getnames <- function(x) sapply(x, attr, "name") 


# see if all the items in a list have the same number of rows 
.equalrows <- function(res) {
  rows <- sapply(res, nrow)
  length(unique(rows)) == 1
}


# see if all the items in a list have the same number of columns
.equalcols <- function(res) {
  cols <- sapply(res, ncol)
  length(unique(cols)) == 1
}


# result directory exists? 
.checkResultDir <- function(resultDir) !is.na(resultDir) & dir.exists(resultDir)
