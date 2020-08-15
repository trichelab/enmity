# BiocManager::install("trichelab/enmity")
library(enmity) 
library(HDF5Array)
library(MultiAssayExperiment) 

if (!exists("meth_hdf5")) meth_hdf5 <-loadHDF5SummarizedExperiment("scNMT_meth")
if (!exists("acc_hdf5")) acc_hdf5 <- loadHDF5SummarizedExperiment("scNMT_acc")
if (!exists("rna_hdf5")) rna_hdf5 <- loadHDF5SummarizedExperiment("scNMT_rna")

assayList <- list(meth = meth_hdf5, acc = acc_hdf5, rna = rna_hdf5)
ExpList <- ExperimentList(assayList)

commonColDat <- Reduce(intersect, 
                       lapply(ExpList, function(x) names(colData(x))))
colDat <- DataFrame(do.call(cbind, 
                            lapply(ExpList, 
                                   function(x) colData(x)[[commonColDat]])))

# the idiot's version of gapped LCS; should make this a function & document it
LCS <- function(...) {
  paste(Reduce(intersect, strsplit(unlist(...), "")), collapse="")
}

# this gets tricky because there are 8 misses
rownames(colDat) <- apply(colDat, 1, LCS)

  # sample mappings: 
  maplist <- apply(colDat, 2, 
                   function(x) DataFrame(primary=names(x), colname=x))
  sampMap <- listToMap(maplist)

  # disambiguate the column names for the main MAE:
  names(colDat) <- paste(commonColDat, names(colDat), sep=".") 

  # assemble the MultiAssayExperiment:
  scNMT <- MultiAssayExperiment(experiments = ExpList,
                                colData = colDat, 
                                sampleMap = sampMap)

  # save that 
  saveRDS(scNMT, file="scNMT.rds") 

  # load it 
  scNMT0 <- readRDS("scNMT.rds") 

  # test it 
  stopifnot(identical(scNMT, scNMT0))

  # move it around
  tmpbase <- basename(scratch)
  newtmp <- sub(tmpbase, reverse(tmpbase), scratch)
  dir.create(newtmp)
  newAssayList <- list() 
  for (hdf5dir in c("scNMT_meth", "scNMT_acc", "scNMT_rna")) {
    asy <- sub("scNMT_", "", hdf5dir)
    asypath <- file.path(newtmp, hdf5dir)
    dir.create(asypath)
    for (hdf5file in c("assays.h5", "se.rds")) {
      file.copy(from=file.path(scratch, hdf5dir, hdf5file),
                to=file.path(newtmp, hdf5dir, hdf5file))
    }
    newAssayList[[asy]] <- loadHDF5SummarizedExperiment(asypath)
  }

  newExpList <- ExperimentList(newAssayList)
  scNMT1 <- MultiAssayExperiment(experiments = newExpList,
                                 colData = colDat, 
                                 sampleMap = sampMap)
  stopifnot(!identical(scNMT0, scNMT1)) # better not be!
   
}

# return to the original location (prior to using a temporary directory)
setwd(oldwd)
