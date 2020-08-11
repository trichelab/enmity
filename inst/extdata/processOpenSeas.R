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
            stop("object does not contain either 'Beta' amongst assay slots")
          })

# for now: 
strategy <- "parallel"
backend <- "inmemory"
# can do EITHER HDF5, OR parallel, but not both, in the large (120 runs) 

tsvs <- list.files(patt="^GS.*_CpG\\-met_processed\\.tsv(\\.gz)?$")
tsvs <- unique(sub("\\.gz", "", tsvs))
if (testing) {

  tsvs <- tsvs[1:3]

  message("Testing serial... ")
  se <- mergeScNMT(tsvs, saveGR=FALSE, BPPARAM=SerialParam())
  stopifnot(!all(is.na(getBeta(se))))
  message("OK.\n\n")
  
  message("Testing multicore... ")
  inmemSE <- mergeScNMT(tsvs, saveGR=FALSE, BPPARAM=MulticoreParam())
  stopifnot(!all(is.na(getBeta(inmemSE))))
  message("OK.\n\n")

  message("Testing multicore HDF5... (Note: this does not scale up properly.)")
  hdf5SE <- mergeScNMT(tsvs, saveGR=FALSE, HDF5=TRUE, BPPARAM=MulticoreParam())
  stopifnot(!all(is.na(getBeta(hdf5SE))))
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
