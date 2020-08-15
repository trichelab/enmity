# if minfi is loaded, this will already be defined: 
if (!exists("getBeta") | !isGeneric("getBeta")) {
  setGeneric("getBeta", function(object, ...) standardGeneric("getBeta"))
}
setMethod("getBeta", signature(object = "SummarizedExperiment"), 
          function(object, ...) {
            nms <- names(assays(object, withDimnames = FALSE))
            if ("Beta" %in% nms) return(assay(object, "Beta"))
            stop("object does not contain 'Beta' amongst assay slots")
          })

# like getBeta but for accessibility 
if (!exists("getAcc") | !isGeneric("getAcc")) {
  setGeneric("getAcc", function(object, ...) standardGeneric("getAcc"))
}
setMethod("getAcc", signature(object = "SummarizedExperiment"), 
          function(object, ...) {
            nms <- names(assays(object, withDimnames = FALSE))
            if ("Acc" %in% nms) return(assay(object, "Acc"))
            stop("object does not contain 'Acc' amongst assay slots")
          })


# like getAcc but for counts
if (!exists("counts") | !isGeneric("counts")) {
  setGeneric("counts", function(object, ...) standardGeneric("counts"))
}
setMethod("counts", signature(object = "SummarizedExperiment"), 
          function(object, ...) {
            nms <- names(assays(object, withDimnames = FALSE))
            if ("counts" %in% nms) return(assay(object, "counts"))
            stop("object does not contain 'counts' amongst assay slots")
          })
