#' take an scNMT CpG methylation file and read its coordinates into a GRanges 
#' 
#' This function is primarily to Reduce(union, aBunchOfGRs) for aggregation
#' 
#' @param filename  the file 
#' @param verbose   squawk? (TRUE) 
#' @param is0based  is it a UCSC file? (FALSE) 
#' @param which     Optional instance of 'GRanges' or 'IntegerRangesList' (NULL)
#' @param ...       Optional arguments, currently ignored
#' 
#' @return          a GRanges with the tabix-derived coordinates, but no $score
#' 
#' @import Rsamtools
#' @import GenomicRanges
#' 
#' @export
scNMTtoGR <- function(filename, verbose=TRUE, is0based=FALSE, which=NULL, ...) {

  tf <- Rsamtools::TabixFile(filename, ...)
  if (!file.exists(Rsamtools::index(tf))) stop("Your file lacks an index!")
  if (verbose) message("Reading GRanges from ", basename(path(tf)))

  chrs <- Rsamtools::headerTabix(tf)$seqnames
  cols <- Rsamtools::headerTabix(tf)$indexColumns
  what <- as.list(rep("character", max(cols)))
  names(what) <- rep("ignore", length(what))
  if (length(what) < length(cols) | (max(cols) < 3)) {
    what <- list("character", "integer", "integer")
    names(what) <- c("seqnames", "start", "end")
  } else { 
    what[cols[2:3]] <- "integer" # start, end
    names(what)[cols] <- c("seqnames", "start", "end")
  }

  tfdf <- as.data.frame(scan(path(tf), what=what, sep="\t", flush=TRUE))[, cols]
  names(tfdf) <- names(what)
  tfdf$start <- as.integer(tfdf$start) # some GTFs have bogus start/end ranges
  tfdf$end <- as.integer(tfdf$end)     # so we coerce them here and then subset
  tfdf <- subset(tfdf, !is.na(seqnames) & !is.na(start) & !is.na(end))
  if (!is.null(which)) {
    tfdf <- subset(tfdf, seqnames %in% unique(seqnames(which)))
  }
  gr <- makeGRangesFromDataFrame(tfdf, starts.in.df.are.0based=is0based)
  if (!is.null(which)) {
    gr <- subsetByOverlaps(gr, which)
  }
  metadata(gr)$filename <- filename
  return(gr)

}
