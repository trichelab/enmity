#' somewhat of an evolved version of tabixToGR with checkpointing
#' 
#' Note that the parameter `which` behaves as with scanBamParam() or import(). 
#' 
#' @param tfl       a TabixFileList or something that can become one
#' @param BPPARAM   BiocParallel parameters (default is SerialParam())
#' @param which     optional GRanges to restrict reading of coordinates (NULL)
#' @param verbose   be verbose? (TRUE) 
#' 
#' @return          a GRanges
#' 
#' @import BiocGenerics
#' @import GenomicRanges
#' 
#' @export
scNMTFileListToGR <- function(tfl, BPPARAM=SerialParam(), which=NULL, verbose=TRUE) {

  # sanity checking per usual 
  if (!is(tfl, "TabixFileList")) {
    tfl <- try(TabixFileList(lapply(tfl, TabixFile)))
    if (inherits(tfl, "try-error")) stop("Invalid tabix files; cannot proceed.")
  }

  # bit of housekeeping to avoid a disaster 
  for (bpdir in c(bplogdir(BPPARAM), bpresultdir(BPPARAM))) {
    if (!is.na(bpdir) & !dir.exists(bpdir)) {
      dir.create(bpdir)
    }
  }

  # ok now let's get to work 
  if (verbose) message("Reading indices from ", length(tfl), " Tabix files...")
  grs <- bptry(bplapply(X=tfl, FUN=scNMTtoGR, BPPARAM=BPPARAM, which=which))
  if (!all(bpok(grs))) {
    stop("scNMTtoGR() encountered errors for these files:\n  ",
         paste(sapply(tfl, path)[!bpok], collapse = "\n  "))
  }

  # are the results saved?  Load them first, if so.
  resultDir <- bpresultdir(BPPARAM)
  if (!is.na(resultDir)) {
    if (verbose) message("Loading saved results from ", resultDir)
    grs <- loadBpFiles(BPPARAM)
    if (verbose) message("OK.")
  }

  # turn results into a GRangesList
  filenames <- sapply(grs, function(x) metadata(x)$filename)
  grl <- GRangesList(grs)[order(sapply(grs, length), decreasing=TRUE)]
  metadata(grl)$filenames <- filenames
  rm(grs)

  # and reduce to a GRanges
  if (verbose) message("Taking the union of ", length(grl), " indices...")
  res <- Reduce(GenomicRanges::union, grl)
  metadata(res)$filenames <- filenames
  return(res)

}
