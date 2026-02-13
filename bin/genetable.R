#!/usr/bin/env Rscript

# ~~ Code description ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# generate a meth table based on genes and CGI rather then DMRs
# (1) DMG-table should contain info about which genomic regions DMRs
# are present in
# (2) CGI table should contain overlapping DMRs, these DMRs can be pulled from 
# DMR-table to identify genes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(genomation)
  library(readr)
  library(scales)
})

ens <- read.csv("data/ensembl_table.csv.gz")
meth <- read.csv("data/PFOS_MCF10A_DMR.csv.gz")

# Import refseq
file <- file.path("data", "refseq_UCSC_GRCh38.bed")
annotations <- readTranscriptFeatures(file)

# Import CpG-island data
file <- file.path("data", "cpgislands_GRCh38.bed")
cgi <- readGeneric(file, keep.all.metadata = TRUE)

# only use well annotated chromsomes
chrom <- c(1:22, "X", "Y")
seqlevels(cgi, pruning.mode = "coarse") <- seqlevels(cgi)[seqlevels(cgi) %in% paste0("chr", chrom)]
cgi$id <- paste0("cgi:", 1:length(cgi))
annotations$cgi <- cgi

### Calculate DMR overlap with genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

genes <- unique(ens$gene_id)

build_genetable <- function(gene){
  cat(paste("Testing", gene, "| # =", which(genes %in% gene)), "\n")

  gene_info <- ens[ens$gene_id %in% gene, ]
  meth_info <- meth[meth$chr %in% gene_info$chr, ]

  gr_gene <- as(gene_info[, 1:4], "GRanges")
  gr_meth <- as(meth_info[, 1:3], "GRanges")

  rows <- data.frame(findOverlaps(gr_gene, gr_meth))

  d <- meth_info[rows$subjectHits, ]

  if ( any(d$sign) ) {
    tmp <- d[d$sign == TRUE, ]

    # output table
    out <- gene_info

    out$DMG <- TRUE
    out$num_cpg <- sum(tmp$feature == "cpg")
    out$num_tile <- sum(tmp$feature == "tile100")

    out$dmr_in_cgi <- sum(!is.na(tmp$cgi_id))
    out$dmr_in_promoter <- sum(tmp$region == "promoter")
    out$dmr_in_exon <- sum(tmp$region == "exon")
    out$dmr_in_intron <- sum(tmp$region == "intron")
    out$dmr_in_intragenic <- sum(tmp$region == "intragenic")

    out$dmr_id_cgi <- paste(tmp$dmr_id[!is.na(tmp$cgi_id)], collapse = ",")
    out$dmr_id_promoter <- paste(tmp$dmr_id[tmp$region == "promoter"], collapse = ",")
    out$dmr_id_exon <- paste(tmp$dmr_id[tmp$region == "exon"], collapse = ",")
    out$dmr_id_intron <- paste(tmp$dmr_id[tmp$region == "intron"], collapse = ",")
    out$dmr_id_intragenic <- paste(tmp$dmr_id[tmp$region == "intragenic"], collapse = ",")

  }else{

    out <- cbind(
      gene_info,
      data.frame(FALSE, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,NA,NA)
    )

    colnames(out)[ (ncol(gene_info) + 1) : ncol(out) ] <- c(
      "DMG", "num_cpg", "num_tile", 
      "dmr_in_cgi", "dmr_in_promoter", "dmr_in_exon", "dmr_in_intron", "dmr_in_intragenic", 
      "dmr_id_cgi", "dmr_id_promoter", "dmr_id_exon", "dmr_id_intron", "dmr_id_intragenic"
    )
  }

  return(out)
}

l1 <- lapply(genes, build_genetable)
genes <- Reduce(function(x,y) rbind(x,y), l1)


# = save results ============================================================= #

filename <- "PFOS_MCF-10A_DMG.Rds"
saveRDS(genes, file = file.path("data", filename))

cat(paste(
  "\n~~ genetable.R complete ~~~~~~~~~~~~~~~~~~~~~~~~~~",
  "\n > Output\t:", filename,
  "\n > Num. DMGs\t:", comma(nrow(genes)),
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
))