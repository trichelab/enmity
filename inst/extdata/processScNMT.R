# BiocManager::install("trichelab/enmity")
library(enmity) 
library(HDF5Array) # for explicit HDF5 loads/saves

# where the sample data is kept:
enmity_ext <- system.file("extdata", package="enmity")
scratch <- tempdir()
oldwd <- getwd() 
setwd(scratch)

# catalog the extdata files for testing 
meth_tsv <- file.path(enmity_ext, list.files(enmity_ext, patt="CpG.*.gz$"))
acc_tsv <- file.path(enmity_ext, list.files(enmity_ext, patt="GpC.*.gz$"))
rna_tab <- file.path(enmity_ext, list.files(enmity_ext, patt="EB_RNA_counts"))

# first test:
message("Testing scanScNMT...") 
scanScNMT(meth_tsv[1]) # test scanScNMT function

# for serial vs parallel:
library(BiocParallel)
BPMULTI <- MulticoreParam()

message("Testing serial meth merge... ")
se <- mergeScNMT(meth_tsv, saveGR=FALSE)
stopifnot(!all(is.na(getBeta(se))))
message("OK.\n\n")
  
message("Testing multicore meth merge... ")
inmemSE <- mergeScNMT(meth_tsv, saveGR=FALSE, BPPARAM=BPMULTI)
stopifnot(!all(is.na(getBeta(inmemSE))))
message("OK.\n\n")

message("Testing equivalence...") 
stopifnot(identical(se, inmemSE))
message("OK.\n\n")

# not sure why, but this crashes Singularity
if (Sys.getenv("SINGULARITY_NAME") == "") {
  message("Testing HDF5 meth merge... ")
  hdf5SE <- mergeScNMT(meth_tsv, saveGR=FALSE, HDF5=TRUE)
  stopifnot(!all(is.na(getBeta(hdf5SE))))
  message("OK.\n\n")
} else {
  # as an alternative, just use saveHDF5SummarizedExperiment on in-memory SEs 
  hdf5SE <- saveHDF5SummarizedExperiment(inmemSE, dir="scNMT_meth",replace=TRUE)
}

message("Testing serial accessibility merge... ")
acc_se <- mergeScNMT(acc_tsv, saveGR=FALSE, what="acc")
stopifnot(!all(is.na(getAcc(acc_se))))
message("OK.\n\n")
  
message("Testing multicore accessibility merge... ")
acc_inmem <- mergeScNMT(acc_tsv, saveGR=FALSE, what="acc", BPPARAM=BPMULTI)
stopifnot(!all(is.na(getAcc(acc_inmem))))
message("OK.\n\n")

message("Testing equivalence...") 
stopifnot(identical(acc_se, acc_inmem))
message("OK.\n\n")

# this crashes singularity otherwise 
if (Sys.getenv("SINGULARITY_NAME") == "") {
  message("Testing HDF5 accessibility merge... ")
  acc_hdf5 <- mergeScNMT(acc_tsv, saveGR=FALSE, what="acc", HDF5=TRUE)
  stopifnot(!all(is.na(getAcc(acc_hdf5))))
  message("OK.\n\n")
} else { 
  # as an alternative, just use saveHDF5SummarizedExperiment on in-memory SEs 
  acc_hdf5 <- saveHDF5SummarizedExperiment(acc_inmem, dir="scNMT_acc",replace=T)
}

# load RNA -- this is a lot easier than the CpG and GpC data
message("Testing RNA load...") # wrap this anyway, though
rna <- read.table(rna_tab, header=TRUE, row=1)
head(rna)

# the counts are against ENSEMBL 87, so let's get their coordinates: 
library(AnnotationHub)
ah <- AnnotationHub()
mmq <- query(ah, "Mus_musculus.GRCm38.87.gtf")
stopifnot(length(mmq) == 1)
mm87 <- mmq[[1]]
# snapshotDate(): 2020-04-27
# downloading 1 resources
# retrieving 1 resource
#  |======================================================================| 100%
# loading from cache
# Importing File into R ..

# get the subset of genes represented in EB_P2D12 and EB_P2F12
rr <- subset(mm87, gene_id %in% rownames(rna) & type == "transcript")
columns <- c("gene_id",
             "gene_name",
             "transcript_id",
             "transcript_name",
             "transcript_biotype")
mcols(rr) <- mcols(rr)[, columns]
rr <- split(rr, rr$gene_id)[rownames(rna)]
cdata <- DataFrame(sample=colnames(rna))
library(SummarizedExperiment)
rna_se <- SummarizedExperiment(assays=list(counts=rna), 
                               rowRanges=rr,
                               colData=cdata)
show(rna_se)

# save it as an HDF5-backed SummarizedExperiment
rna_hdf5 <- saveHDF5SummarizedExperiment(rna_se, dir="scNMT_rna")

# wrap up a MultiAssayExperiment
library(MultiAssayExperiment) 
assayList <- list(meth = hdf5SE, 
                  acc = acc_hdf5,
                  rna = rna_hdf5)
ExpList <- ExperimentList(assayList)
MAE <- MultiAssayExperiment(experiments = ExpList,
                            colData = colDat, 
                            sampleMap = sampMap)

# save that 
saveRDS(MAE, file="MAE.rds") 

# load it 
MAE1 <- readRDS("MAE.rds") 

# test it 
stopifnot(identical(MAE, MAE1))

# optional: relocate an HDF5-backed version to test Marcel's debugging?
MAEtest <- FALSE
if (MAEtest) {
  tmpbase <- basename(scratch)
  newtmp <- sub(tmpbase, reverse(tmpbase), scratch)
  dir.create(newtmp)
  newAssayList <- list() 
  for (hdf5dir %in% c("scNMT_meth", "scNMT_acc", "scNMT_rna")) {
    asy <- sub("scNMT_", "", hdf5dir)
    asypath <- file.path(newtmp, hdf5dir)
    dir.create(asypath)
    for (hdf5file %in% c("assays.h5", "se.rds")) {
      file.copy(from=file.path(scratch, hdf5dir, hdf5file),
                to=file.path(newtmp, hdf5dir, hdf5file))
    }
    newAssayList[[asy]] <- loadHDF5SummarizedExperiment(asypath)
  }
  newExpList <- ExperimentList(newAssayList)

}

# return to the original location (prior to using a temporary directory)
setwd(oldwd)
