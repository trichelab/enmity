#' merge a bunch of scNMT raw (met|acc) files (chrom, start, %meth/%acc format)
#'
#' NOTE: see bsseq:::.constructCounts for some background on doing this in HDF5
#' 
#' @param tsvs    scNMT raw CpG (methylation) or GpC (accessibility) files 
#' @param gen     what genome they correspond to (default is GRCm38)
#' @param loci    a pre-merged GRanges with the union of all tsvs' ranges (NULL)
#' @param saveGR  save the intermediate `loci` object, if loci=NULL? (TRUE) 
#' @param saveSE  save the intermediate `se` object before filling it? (FALSE) 
#' @param HDF5    use an HDF5-backed SummarizedExperiment for loading? (FALSE)
#' @param what    methylation ("meth") or accessibility ("acc") data? ("meth")
#' @param dir     if using HDF5 and/or saving loci, path to use ("scNMT_"what)
#' @param seqinf  add seqInfo? (TRUE; set FALSE if running behind a firewall) 
#' @param verbose be extra verbose? (FALSE) 
#' @param which   optional GRanges to restrict reading of coordinates (NULL)
#' @param BPPARAM BiocParallelParam object (default BiocParallel::SerialParam())
#'
#' @return        a SummarizedExperiment, perhaps HDF5-backed, with merged data
#' 
#' @import SummarizedExperiment
#' @import GenomeInfoDb
#' @import BiocParallel 
#' @import DelayedArray
#' @import HDF5Array
#' @import S4Vectors
#' @import bsseq
#' 
#' @export
mergeScNMT <- function(tsvs, gen="GRCm38", loci=NULL, saveGR=TRUE, saveSE=FALSE, HDF5=FALSE, what=c("meth", "acc"), dir="scNMT", seqinf=TRUE, verbose=FALSE, which=NULL, BPPARAM=SerialParam()) { 

  what <- match.arg(what) 
  dir <- paste(dir, what, sep="_") 
  asy <- switch(what, "meth"="Beta", "acc"="Acc")

  # can use HDF5 OR parallel, but not both 
  if (!is(BPPARAM, "SerialParam") & HDF5 == TRUE) { 
    message("You may use parallel processing OR HDF5 backing, but not both.")
    message("Switching to serial processing with an HDF5 backend...") 
    BPPARAM <- SerialParam()
  }

  # (try to) ensure TSVs are tabixed by dry-running scanScNMT
  if (is.null(loci)) { 
    if (verbose) message("Cataloging loci...")
    indexed <- bptry(bplapply(X=tsvs, FUN=scanScNMT, dry=TRUE, BPPARAM=BPPARAM))
    if (!all(bpok(indexed))) { 
      # {{{ usually a fixable issue: the scNMT headers on GEO are fugged
      message("If you see Tabix errors, try (at a shell prompt, not in R):")
      message("") 
      message("for i in GSM*.tsv.gz; do")
      message("  j=`basename $i .gz`; ")
      message("  zcat $i | grep -v rate > $j && \\")
      message("  bgzip -f $j && \\")
      message("  tabix -s 1 -b 2 -e 2 $j.gz; ")
      message("done") 
      message("") 
      # }}}
    }
  }

  # now give them simplified column names
  tsvs <- sub(".gz", "", fixed=TRUE, tsvs)
  tsvpatt <- switch(what,
                    "acc"="_GpC-acc_processed.tsv", 
                    "meth"="_CpG-met_processed.tsv")
  tsvnames <- sub(tsvpatt, "", fixed=TRUE, basename(tsvs))
  tsvnames <- sub("^GSM[0123456789]+_", "", tsvnames) 
  names(tsvs) <- tsvnames
  tsvgzs <- paste0(tsvs, ".gz")
  names(tsvgzs) <- tsvnames
  
  # create a union'ed rowRanges for the se if not provided
  # refactor this?
  if (is.null(loci)) {
    message("Constructing merged locus index `loci`.")
    # Pete suggests to try this with bpiterate instead?
    loci <- scNMTFileListToGR(tsvgzs, verbose=verbose, 
                              BPPARAM=BPPARAM, which=which)
    names(loci) <- as.character(loci) 
    metadata(loci)$files <- tsvgzs
    if (seqinf) {
      message("Adding seqinfo...")
      loci <- .addSeqinfo(loci, gen=gen)
    }
    message("Finished merging into `loci`.")
    if (saveGR) {
      if (!dir.exists(dir)) dir.create(dir)
      lociRDS <- file.path(dir, "loci.rds")
      message("Saving `loci` to ", lociRDS, " ... ", appendLF=FALSE)
      saveRDS(loci, file=lociRDS)
      message("OK.") 
    }
  } else {
    if (!is(loci, "GRanges")) stop("`loci` must be a GRanges, if not NULL.")
    if (is.null(names(loci))) names(loci) <- as.character(loci)
  }

  # create a SummarizedExperiment to hold the data
  message("Creating a SummarizedExperiment for metadata.")
  cdata <- DataFrame(sample=tsvnames, file=tsvgzs)
  se <- SummarizedExperiment(rowRanges=loci, colData=cdata)
  stopifnot(all(file.exists(se$file)))
  colnames(se) <- cdata$sample
  genome(se) <- gen

  # save empty `se`, if requested
  if (saveSE) {
    if (!dir.exists(dir)) dir.create(dir)
    emptySeRDS <- file.path(dir, "empty_se.rds")
    message("Saving empty `se` to ", emptySeRDS , " ... ", appendLF=FALSE)
    saveRDS(se, file=emptySeRDS)
    message("OK.") 
  }

  # NOTE: a lot of the following is straight-up ripped off from bsseq!
  # If we had both methylated read counts and total read counts for scNMT data,
  # we could just use bsseq for the backend and elide most of this! (probably)
  ans_nrow <- length(loci)
  ans_ncol <- length(tsvgzs) 
  ans_dim <- c(ans_nrow, ans_ncol)
  # NOTE: should we use h5writeDimnames to record the dimnames?
  
  dat_type <- "double"
  mat_type <- paste(ifelse(HDF5, "HDF5", "in-memory"), dat_type, "matrix")
  # allocate a big enough matrix (may switch to tiledb instead of HDF5)
  message("Allocating a ", ans_nrow, " x ", ans_ncol, " ", mat_type, ".")
  grid <- RegularArrayGrid(refdim = ans_dim, spacings = c(ans_nrow, 1L))
  DelayedArray:::set_verbose_block_processing(TRUE)
  message("Creating a DelayedArray for the data.")
  resultDir <- bpresultdir(BPPARAM) # if saving

  if (HDF5) {

    message("Warning: HDF5 writes often crash under Singularity.") 

    if (!dir.exists(dir)) dir.create(dir) 
    h5_path <- file.path(dir, "assays.h5")
    if (file.exists(h5_path)) unlink(h5_path)
    asy_sink <- HDF5RealizationSink(dim = ans_dim,
                                    type = "double",
                                    filepath = h5_path,
                                    name = asy)
    on.exit(close(asy_sink), add = TRUE)

    sink_lock <- ipcid()
    on.exit(ipcremove(sink_lock), add = TRUE)

  } else {

    asy_sink <- NULL
    sink_lock <- NULL

  }

  # read in the asy values from scNMT files
  if (verbose) message("Reading in the assay data...")
  asy_dat <- bptry(bplapply(X = seq_along(grid),
                            FUN = .updateScNMT, 
                            files = tsvgzs,
                            loci = rowRanges(se),
                            grid = grid,
                            asy_sink = asy_sink,
                            sink_lock = sink_lock,
                            result_dir = resultDir,
                            gen = gen, 
                            BPPARAM = BPPARAM))

  # checkpoint: 
  if (!all(bpok(asy_dat))) { 
    stop(".updateScNMT() encountered errors for these files:\n  ",
         paste(files[!bpok], collapse = "\n  "))
  }

  # are the results saved?  Load them first, if so.
  if (!is.na(resultDir)) {
    if (verbose) message("Loading saved results from ", resultDir)
    asy_dat <- .loadBpFiles(BPPARAM)
    if (verbose) browser() # check naming
    if (verbose) message("OK.")
  }

  
  # write 'em  
  if (HDF5) {

    asy_dat <- as(asy_sink, "DelayedArray")
    stopifnot(identical(dim(asy_dat), dim(se)))
    assay(se, asy, withDimnames=FALSE) <- asy_dat
    x <- se 
    x@assays <- HDF5Array:::.shorten_assay2h5_links(x@assays)
    saveRDS(x, file = file.path(dir, "se.rds"))

  } else { 
    
    asy_dat <- Reduce(cbind, asy_dat)
    stopifnot(identical(attr(asy_dat, "dim"), ans_dim))
    rownames(asy_dat) <- names(loci)
    colnames(asy_dat) <- colnames(se)
    assay(se, asy) <- asy_dat

  } 

  # done
  return(se)

} 


# utility fn, stolen from bsseq, more or less
.updateScNMT <- function(i, files, loci, grid, asy_sink, sink_lock, gen, result_dir=NULL, which=NULL) {

  if (!is.null(result_dir)) stopifnot(!is.null(names(files))) # need names
  message("[.updateScNMT] Extracting scores for ", names(files)[i])
  message("               from ", files[i])
  gr <- scanScNMT(files[i], gen = gen, which=which)
  ol <- findOverlaps(gr, loci) # does this need to be `equal`?!
  asy_dat <- matrix(rep(NA_real_, length(loci)), ncol = 1)
  asy_dat[subjectHits(ol)] <- score(gr[queryHits(ol)])
  attr(asy_dat, "name") <- names(files)[i]  
  attr(asy_dat, "file") <- files[i]
  if (is.null(asy_sink)) return(asy_dat) # in-memory, perhaps saved

  message("[.updateScNMT] Locking and writing DelayedArray...")
  viewport <- grid[[i]] # HDF5 or similar backend
  ipclock(sink_lock) # respect locking 
  write_block(x = asy_sink, viewport = viewport, block = asy_dat)
  ipcunlock(sink_lock) # respect locking
  message("               Written and unlocked.")
  return(NULL)

}
