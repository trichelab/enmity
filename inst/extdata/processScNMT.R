# BiocManager::install("trichelab/enmity")
library(enmity) 

# we don't typically test on HPC 
host <- system2("hostname", stdout = T)
if (!exists("testing")) testing <- !grepl("triche|node", host)

if (testing) { 
  # arrives from GEO in gzipped format; need to bgzip it instead
  tsv <- "GSM2936019_EBcontrol_P1D12_CpG-met_processed.tsv" 
  scanScNMT(tsv) # test scanScNMT function
}

# if minfi is loaded, this will already be defined: 
if (!exists("getBeta") | !isGeneric("getBeta")) {
  setGeneric("getBeta", function(object, ...) standardGeneric("getBeta"))
}

# could also specialize this for bsseq objects and the like
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

# could also specialize this for bsseq objects and the like
setMethod("getAcc", signature(object = "SummarizedExperiment"), 
          function(object, ...) {
            nms <- names(assays(object, withDimnames = FALSE))
            if ("Acc" %in% nms) return(assay(object, "Acc"))
            stop("object does not contain 'Acc' amongst assay slots")
          })


# for now: 
strategy <- "parallel"
backend <- "inmemory"
# can do EITHER HDF5, OR parallel, but not both, in the large (120 runs) 

meth_tsvs <- list.files(patt="^GS.*_CpG\\-met_processed\\.tsv(\\.gz)?$")
meth_tsvs <- unique(sub("\\.gz", "", meth_tsvs))

acc_tsvs <- list.files(patt="^GS.*_GpC\\-acc_processed\\.tsv(\\.gz)?$")
meth_tsvs <- unique(sub("\\.gz", "", acc_tsvs))

if (testing) {

  tsvs <- meth_tsvs[1:2]

  message("Testing serial... ")
  se <- mergeScNMT(tsvs, saveGR=FALSE, BPPARAM=SerialParam())
  stopifnot(!all(is.na(getBeta(se))))
  message("OK.\n\n")
  
  message("Testing multicore... ")
  inmemSE <- mergeScNMT(tsvs, saveGR=FALSE, BPPARAM=MulticoreParam())
  stopifnot(!all(is.na(getBeta(inmemSE))))
  message("OK.\n\n")

  message("Testing HDF5... ")
  hdf5SE <- mergeScNMT(tsvs, saveGR=FALSE, HDF5=TRUE)
  stopifnot(!all(is.na(getBeta(hdf5SE))))
  message("OK.\n\n")

  accs <- acc_tsvs[1:2]

  message("Testing serial... ")
  acc_se <- mergeScNMT(accs, saveGR=FALSE, what="acc", BPPARAM=SerialParam())
  stopifnot(!all(is.na(getAcc(acc_se))))
  message("OK.\n\n")
  
  message("Testing multicore... ")
  acc_inmem <- mergeScNMT(accs,saveGR=FALSE,what="acc",BPPARAM=MulticoreParam())
  stopifnot(!all(is.na(getAcc(acc_inmem))))
  message("OK.\n\n")

  message("Testing HDF5... ")
  acc_hdf5 <- mergeScNMT(accs, saveGR=FALSE, what="acc", HDF5=TRUE)
  stopifnot(!all(is.na(getAcc(acc_hdf5))))
  message("OK.\n\n")

} else { 

  library(HDF5Array)
  loci <- NULL
  lociRDS <- file.path("scNMT_meth", "loci.rds")
  if (file.exists(lociRDS)) {
    message("Loading pre-saved loci...")
    loci <- readRDS(lociRDS)
  }

  # for stability purposes 
  HDF5 <- FALSE 
  if (backend == "HDF5") {
    HDF5 <- TRUE
    strategy <- "serial"
  }
  BPPARAM <- switch(strategy,
                    `serial`=SerialParam(), 
                    `parallel`=MulticoreParam())
  se <- mergeScNMT(tsvs, loci=loci, HDF5=HDF5, BPPARAM=BPPARAM)
  # Passively saved as HDF5 if HDF5 == TRUE 
  
  if (backend == "inmemory") {
    message("Saving merged scNMT methylation data as HDF5... ", appendLF=FALSE)
    se <- saveHDF5SummarizedExperiment(se, dir="scNMT_meth", replace=TRUE)
    message("done.")
  }

}

show(se)
