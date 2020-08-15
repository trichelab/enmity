# BiocManager::install("trichelab/enmity")
library(enmity) 
enmity_ext <- system.file("extdata", package="enmity")

# catalog the extdata files for testing 
meth_tsv <- file.path(enmity_ext, list.files(enmity_ext, patt="CpG.*.gz$"))
acc_tsv <- file.path(enmity_ext, list.files(enmity_ext, patt="GpC.*.gz$"))
rna_tab <- file.path(enmity_ext, list.files(enmity_ext, patt="EB_RNA_counts"))

# first test:
message("Testing scanScNMT...") 
scanScNMT(meth_tsv[1]) # test scanScNMT function

# for serial vs parallel:
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

message("Testing HDF5 meth merge... ")
hdf5SE <- mergeScNMT(meth_tsv, saveGR=FALSE, HDF5=TRUE)
stopifnot(!all(is.na(getBeta(hdf5SE))))
message("OK.\n\n")

message("Testing serial accessibility merge... ")
acc_se <- mergeScNMT(acc_tsv, saveGR=FALSE, what="acc")
stopifnot(!all(is.na(getAcc(acc_se))))
message("OK.\n\n")
  
message("Testing multicore... ")
acc_inmem <- mergeScNMT(acc_tsv, saveGR=FALSE, what="acc", BPPARAM=BPMULTI)
stopifnot(!all(is.na(getAcc(acc_inmem))))
message("OK.\n\n")

message("Testing equivalence...") 
stopifnot(identical(acc_se, acc_inmem))
message("OK.\n\n")
  
message("Testing HDF5... ")
acc_hdf5 <- mergeScNMT(acc_tsv, saveGR=FALSE, what="acc", HDF5=TRUE)
stopifnot(!all(is.na(getAcc(acc_hdf5))))
message("OK.\n\n")

# load RNA -- this is a lot easier than the CpG and GpC data
message("Testing RNA...")
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
cdata <- DataFrame(cell=colnames(rna))
rna_se <- SummarizedExperiment(assays=list(counts=rna), 
                               rowRanges=rr,
                               colData=cdata)
show(rna_se)

# wrap up a MultiAssayExperiment
library(MultiAssayExperiment) 
