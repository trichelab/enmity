#' @export
if (!exists("getBeta") | !isGeneric("getBeta")) {
  setGeneric("getBeta", function(object, ...) standardGeneric("getBeta"))
}


#' @export
setMethod("getBeta", signature(object = "SummarizedExperiment"), 
          function(object, ...) {
            nms <- names(assays(object, withDimnames = FALSE))
            if ("Beta" %in% nms) return(assay(object, "Beta"))
            stop("object does not contain 'Beta' amongst its named assays")
          })


#' @export
if (!exists("getAcc") | !isGeneric("getAcc")) {
  setGeneric("getAcc", function(object, ...) standardGeneric("getAcc"))
}


#' @export
setMethod("getAcc", signature(object = "SummarizedExperiment"), 
          function(object, ...) {
            nms <- names(assays(object, withDimnames = FALSE))
            if ("Acc" %in% nms) return(assay(object, "Acc"))
            stop("object does not contain 'Acc' amongst its named assays")
          })


#' @export
if (!exists("counts") | !isGeneric("counts")) {
  setGeneric("counts", function(object, ...) standardGeneric("counts"))
}


#' @export
setMethod("counts", signature(object = "SummarizedExperiment"), 
          function(object, ...) {
            nms <- names(assays(object, withDimnames = FALSE))
            if ("counts" %in% nms) return(assay(object, "counts"))
            stop("object does not contain 'counts' amongst its named assays")
          })
