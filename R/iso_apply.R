#' sensible row- or column- wise operations against a DelayedMatrix
#' 
#' For good performance, entire chunks should be loaded; consequently, 
#' if MARGIN == 1 (rows), `by` should be a multiple of chunkdim(X)[1], and
#' if MARGIN == 2 (columns), `by` should be a multiple of chunkdim(X)[2]. 
#' (There is a check included for this, but the principle is worth noting.)
#' 
#' Remember to setRealizationBackend("HDF5Array") or TileDBArray or whatever.
#' 
#' @param X       a DelayedArray (e.g. a DelayedMatrix backed by HDF5Array)
#' @param MARGIN  apply across rows (1) or columns (2), as in base::apply
#' @param FUN     the function to apply to each row or column (see above)
#' @param by      the number of rows or columns to load at one time
#'
#' @return        a new DelayedArray (presumably of the same dimensions as X)
#' 
#' @import        DelayedArray
#' 
#' @export
iso_apply <- function(X, MARGIN, FUN, by=NULL)
{   
    FUN <- match.fun(FUN)
    X_dim <- dim(X)
    stopifnot(length(X_dim) == 2L, isSingleNumber(MARGIN), MARGIN %in% 1:2)
    sink <- RealizationSink(X_dim, dimnames=dimnames(X), type=type(X))

    # check for a sane value of `by` 
    if (!is.null(by)) {
        chunksize <- chunkdim(X)[MARGIN]
        if ((by %% chunksize) > 0) {
            message("Warning: `by` (",by,") is not a multiple of ",
                    "chunkdim(X)[", MARGIN, "] (i.e., ", chunksize, ").")
        }
    }

    ## Use row/col[Auto]Grid() to create a grid of blocks
    ## where the blocks are made of full rows/columns.
    newerDelayedArray <- (packageVersion("DelayedArray") >= 0.15)
    rowFn <- ifelse(newerDelayedArray, rowAutoGrid, rowGrid)
    colFn <- ifelse(newerDelayedArray, colAutoGrid, colGrid)
    X_grid <- switch(MARGIN, `1`=rowFn(X, nrow=by), `2`=colFn(X, ncol=by))
    
    nblock <- length(X_grid)
    for (b in seq_len(nblock)) {
        X_block <- read_block(X, X_grid[[b]])
        Y_block <- apply(X_block, MARGIN, FUN)
        write_block(sink, X_grid[[b]], Y_block)
    }
    as(sink, "DelayedArray")
}
