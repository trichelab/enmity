library(enmity)

# GSE109262
countmats <- list.files(patt="GSE109262.*counts.txt.gz")
countses <- lapply(countmats, read.table, head=TRUE, row=1)
stopifnot(identical(rownames(countses[[1]]), rownames(countses[[2]])))
rna <- do.call(cbind, countses)

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
# class: RangedSummarizedExperiment 
# dim: 22084 128 
# metadata(0):
# assays(1): counts
# rownames(22084): ENSMUSG00000051951 ENSMUSG00000025900 ...
#   ENSMUSG00000096730 ENSMUSG00000095742
# rowData names(0):
# colnames(128): EB_P1A01 EB_P1D01 ... ESC_H09 ESC_H10
# colData names(1): sample

# save it as an HDF5-backed SummarizedExperiment
rna_hdf5 <- saveHDF5SummarizedExperiment(rna_se, dir="scNMT_rna", replace=TRUE)

