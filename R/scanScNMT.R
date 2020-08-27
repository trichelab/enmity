#' import an scNMT raw methylation file (chrom, start, %meth)
#'
#' @param tsv     a tsv or tsv.gz filename (will figure out what to do with it)
#' @param dry     dry run? (just check that the file is tabixed?) (FALSE)
#' @param gen     what genome the file corresponds to (default is GRCm38)
#' @param seqinf  add seqInfo? (leave FALSE if running behind a firewall) 
#' @param verbose be verbose? (FALSE) 
#' @param which   optional GRanges to restrict reading of coordinates (NULL)
#' @param ...     additional arguments to pass on to rtracklayer::import
#'
#' @return        a GRanges with appropriate seqinfo and scores turned to betas
#' 
#' @import Rsamtools
#' @import rtracklayer
#' @import GenomeInfoDb
#' @import GenomicRanges
#' 
#' @export
scanScNMT <- function(tsv, dry=FALSE, gen="GRCm38", seqinf=FALSE, verbose=FALSE, which=NULL, ...) {

  if (verbose) message("Scanning ", tsv, "...")
  
  # tidy up the input filename 
  tsv <- sub("(\\.gz)+$", "", tsv)
  if (!grepl("(CpG|GpC)\\-(met|acc)_processed.tsv", tsv)) {
    warning("This does not look like an scNMT CpG or GpC file. It may fail.")
  }

  # ensure it is bgzipped 
  tsvgz <- paste0(tsv, ".gz")
  if (verbose) message("Checking for bgzipped ", tsvgz, "...")
  if (!file.exists(tsvgz)) {
    message("bgzipping...")
    stopifnot(file.exists(tsv))
    tsvgz <- bgzip(tsv, dest=tsvgz)
  } 
  
  # ensure it is indexed 
  tsvtbi <- paste0(tsvgz, ".tbi")
  if (verbose) message("Checking for index ", tsvtbi, "...")
  if (!file.exists(tsvtbi) ) {
    tsvtbi <- try(.indexScNMT(tsvgz, verbose=FALSE), silent=TRUE)
    if (inherits(tsvtbi, "try-error")) {
      # gunzip and then bgzip the tsv file
      if (!file.exists(tsv)) tsv <- gunzip(tsvgz)
      tsvgz <- bgzip(tsv, dest=tsvgz, overwrite=TRUE)
      tsvtbi <- .indexScNMT(tsvgz)
    }
  }

  if (dry) {
    message("Checking index for ", tsvgz, "... ", appendLF=FALSE) 
    stopifnot(is(TabixFile(tsvgz), "TabixFile")) 
    message("OK.")
  } else { 
    gr <- import(TabixFile(tsvgz), which=which)
    names(mcols(gr))[1] <- "score"
    maxscore <- max(gr$score) # typically 0-100 in raw
    if (maxscore > 1) gr$score <- gr$score / maxscore
    names(gr) <- as.character(gr) # for merging
    if (seqinf) gr <- .addSeqinfo(gr, gen=gen)
    return(gr) 
  }

}


# utility fn
.indexScNMT <- function(tsvgz, verbose=TRUE) {
  if (verbose) message("Indexing ", tsvgz, "... ", appendLF=FALSE)
  res <- indexTabix(tsvgz, seq=1, start=2, end=2, skip=-1L)
  if (verbose) message("OK.")
  return(res)
}


# utility fn -- refactor to use h5testR::getSeqinfo
.addSeqinfo <- function(gr, gen=NULL) { 

  if (is.null(gen)) gen <- unique(genome(gr))
  if (is.na(gen) | is.null(gen)) stop("No genome detected!")
  columns <- c("SequenceName", "SequenceLength", "circular")
  seqinf <- with(getChromInfoFromNCBI(gen)[, columns], 
                 Seqinfo(SequenceName, SequenceLength, circular, gen))
  seqinfo(gr) <- seqinf[seqlevels(gr)]
  return(gr)

}
