#!/usr/bin/env Rscript

# ~~ Code description ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Use only DMRs w/ pval < 0.05 & abs(meth.diff) > 15
# Overlap w/ genomic regions
# import all differentally methylated tables from the folder data
# merge into one file
# DMR-table should contain overlapping gene IDs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

suppressPackageStartupMessages({
  library(biomaRt)
  library(GenomicRanges)
  library(genomation)
  library(readr)
  library(scales)
})


ens <- read_csv("data/ensembl_table.csv.gz")
meth <- read_csv("data/diffmeth.csv.gz")

# Add significance column
meth$sign <- FALSE
rows <- meth$qvalue <= 0.05 & abs(meth$meth.diff) >= 15
meth$sign[rows] <- TRUE

# Reduce computation time
meth <- meth[meth$sign == TRUE, ]


# Make column w/ overlapping genes per DMR, multiple genes are comma-separated

# Make a range of all ensembl genes
gr_ens <- as(ens[, 1:4], "GRanges")

# Overlap with range of all DMRs
meth$gene_id <- apply(meth, 1, function(x) {
  #show(x[9])
  gr_meth <- as( as.data.frame(t(x[1:3])), "GRanges")
  rows <- data.frame(findOverlaps(gr_ens, gr_meth))
  #show(nrow(rows))

  if (nrow(rows) == 0) out <- NA
  if (nrow(rows) > 0) out <- paste(ens$gene_id[rows$queryHits], collapse = ",")

  return(out)
})


### CpG-island positions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Import CpG-island data
file <- file.path("data", "cpgislands_GRCh38.bed")
cgi <- readGeneric(file, keep.all.metadata = TRUE)

# only use well annotated chromsomes
chrom <- c(1:22, "X", "Y")
seqlevels(cgi, pruning.mode = "coarse") <- seqlevels(cgi)[seqlevels(cgi) %in% paste0("chr", chrom)]
cgi$id <- paste0("cgi:", 1:length(cgi))

# add cpgisland data, include as ID
gr_meth <- as(meth[, 1:3], "GRanges")
rows <- data.frame(findOverlaps(cgi, gr_meth))
meth$cgi_id <- NA
meth$cgi_id[rows$subjectHits] <- cgi$id[rows$queryHits]


### Genomic features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import refseq
file <- file.path("data", "refseq_UCSC_GRCh38.bed")
annotations <- readTranscriptFeatures(file)

# intergenic
meth$region <- ifelse(is.na(meth$gene_id), "intergenic", "intragenic")

# introns
rows <- data.frame(findOverlaps(annotations$intron, gr_meth))
meth$region[rows$subjectHits] <- "intron"

# exon
rows <- data.frame(findOverlaps(annotations$exon, gr_meth))
meth$region[rows$subjectHits] <- "exon"

# promoter
rows <- data.frame(findOverlaps(annotations$promoter, gr_meth))
meth$region[rows$subjectHits] <- "promoter"

# FIX, if no gene id set to intergenic (between genes)
meth$region <- ifelse(is.na(meth$gene_id), "intergenic", meth$region)

filename <- paste0("PFOS_MCF10A_DMR.csv.gz")
write_csv(meth, file = file.path("data", filename))

cat(paste(
  "\n~~ ensembl.R complete ~~~~~~~~~~~~~~~~~~~~~~~~~~",
  "\n > Output\t:", filename,
  "\n > Num. DMRs\t:", comma(nrow(meth)),
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
))