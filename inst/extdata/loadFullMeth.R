# BiocManager::install("trichelab/enmity")
library(enmity) 
library(HDF5Array) # for HDF5 loads/saves
library(BiocParallel) # for multiprocessing 
BPMULTI <- MulticoreParam()

# catalog the extdata files for testing 
meth_tsv <- list.files(patt="CpG.*.gz$")
meth_se <- mergeScNMT(meth_tsv, saveGR=FALSE, BPPARAM=BPMULTI)
stopifnot(!all(is.na(getBeta(inmemSE))))
meth_hdf5 <- saveHDF5SummarizedExperiment(meth_se, dir="scNMT_meth", replace=T)
