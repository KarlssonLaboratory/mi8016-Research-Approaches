#!/usr/bin/env Rscript

# ~~ Code description ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Gene ontology enrichment analysis
# Divided between subgroups: DMRs in (1) promoter, (2) exon, (3) intron,
# (4) CGI
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# Import datasets
genes <- read.csv("data/PFOS_MCF10A_DMG.csv.gz")
meth <- read.csv("data/PFOS_MCF10A_DMR.csv.gz")
ens <- read.csv("data/ensembl_table.csv.gz")

# Extract genes for each genomic feature: promoter, exon, CpG-island
geneset <- list()

cgi <- na.omit(subset(meth, !is.na(cgi_id))$gene_id)
cgi <- unlist(sapply(cgi, function(x) strsplit(x, ",")[[1]]))
cgi <- unique(unname(cgi))

geneset$promoter$gene_id <- na.omit(subset(genes, dmr_in_promoter > 0)$gene_id)
geneset$promoter$entrez_id <- na.omit(subset(genes, dmr_in_promoter > 0)$entrez_id)
geneset$exon$gene_id <- na.omit(subset(genes, dmr_in_exon > 0)$gene_id)
geneset$exon$entrez_id <- na.omit(subset(genes, dmr_in_exon > 0)$entrez_id)
geneset$cgi$gene_id <- cgi
geneset$cgi$entrez_id <- na.omit(ens[ens$gene_id %in% cgi, "entrez_id"])

go_analysis <- function(genelist, genes){

  cat(" > Running KEGG analysis...\n")

  go_kegg <- clusterProfiler::enrichKEGG(
    gene = genelist$entrez_id,
    organism = "hsa"
  )
  go_kegg <- data.frame(go_kegg)

  cat(" > Running GO analysis...\n")

  go <- clusterProfiler::enrichGO(
    gene = genelist$gene_id,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "ALL",
    universe = genes$gene_id,
    readable = TRUE
  )
  go <- data.frame(go)

  cat(" > Merging GO results...\n")

  # check if any sign GO list was found
  if (nrow(go) > 0 || nrow(go_kegg) > 0) {
    if (nrow(go_kegg) > 0) {
      go_kegg$ONTOLOGY <- "KEGG"

      # translate entrez_id to gene name
      go_kegg$geneID <- sapply(go_kegg$geneID, function(gene){
        genes <- strsplit(gene, "/")[[1]]
        genes <- subset(ens, entrez_id %in% genes)$gene_name
        paste(genes, collapse = "/")
      })

      go <- rbind(go, go_kegg[, 3:ncol(go_kegg)])
    }
  }

  cat(" > Sorting GO results...\n")

  # sort on qvalue & split the gene IDs
  go <- go[order(go$qvalue), ]
  go$ONTOLOGY <- factor(go$ONTOLOGY, levels = c("BP", "CC", "MF", "KEGG"))
  go$geneID <- sapply(go$geneID, function(gene) {
    strsplit(gene, "/")[[1]]
  })
  go$gene_id <- sapply(go$geneID, function(gene) {
    ens[ens$gene_name %in% gene, "gene_id"]
  })

  cat(" > Calculating gene ranks...\n")

  # Rank genes by occurrence in the GO terms
  go_genes <- unique(unlist(go$geneID))
  terms <- unique(go$ONTOLOGY)

  d <- sapply(terms, function(x) {
    tab <- subset(go, ONTOLOGY %in% x)$geneID
    tab <- table(unlist(tab))
    rows <- match(names(tab), go_genes)

    dat <- vector(mode = "numeric", length = length(go_genes))
    names(dat) <- go_genes
    dat[rows] <- tab

    dat
  })
  colnames(d) <- terms

  # sort by freq. of genes in go terms by rank
  ranks <- data.frame(apply(d, 2, rank))
  ranks$tot <- rowSums(ranks)
  rows <- order(-ranks$tot)

  generanks <- data.frame(gene_name = go_genes, d)
  generanks <- generanks[rows,]
  rownames(generanks) <- NULL

  res <- list()
  res$go <- go
  res$generanks <- generanks

  return(res)
}

res <- lapply(geneset, function(genelist) go_analysis(genelist, genes))

filename <- "PFOS_MCF-10A_GO.Rds"
saveRDS(res, file = file.path("data", "PFOS_MCF-10A_GO.Rds"))

cat(paste(
  "\n~~ gene_ontology.R complete ~~~~~~~~~~~~~~~~~~~~~~~~~~",
  "\n > Output\t:", filename,
  "\n > Num. DMGs\t:", comma(length(res)),
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
))