# BiocManager::install("trichelab/enmity")
library(enmity) 
library(HDF5Array) # for HDF5 loads/saves
library(BiocParallel) # for multiprocessing 
BPMULTI <- MulticoreParam()

# catalog the extdata files for testing 
acc_tsv <- list.files(patt="GpC.*.gz$")
acc_se <- mergeScNMT(acc_tsv, saveGR=FALSE, BPPARAM=BPMULTI)
stopifnot(!all(is.na(getAcc(acc_se))))
acc_hdf5 <- saveHDF5SummarizedExperiment(acc_se, dir="scNMT_acc", replace=TRUE)
