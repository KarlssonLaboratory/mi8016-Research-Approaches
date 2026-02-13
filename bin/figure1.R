#!/usr/bin/env Rscript

# ~~ Code description ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Generate figure 1 from paper:
# https://doi.org/10.1016/j.scitotenv.2024.174864
# Figure is composed of five attributes
# 1. Methylation status of DMRs
# 2. Genomic regions of DMRs
# 3. DM CpG-islands
# 4. Genes overlapping DMRs
# 5. Gene types of DMRs overlapping with promoter, exon and/or CGI
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

suppressPackageStartupMessages({
  library(tidyverse)
  library(RColorBrewer)
  library(ggrepel)
  library(scales)
  library(genomation)
  library(GenomicRanges)
  library(cowplot)
})

# Import datasets
meth <- readRDS("data/PFOS_MCF-10A_DMR.Rds")
genes <- readRDS("data/PFOS_MCF-10A_DMG.Rds")
cpg_islands <- readGeneric(
  "data/cpgislands_GRCh38.bed",
  keep.all.metadata = TRUE
)

# Only keep the annotated chromosomes
chrom <- c(1:22, "X")
seqlevels(cpg_islands, pruning.mode = "coarse") <- seqlevels(cpg_islands)[seqlevels(cpg_islands) %in% paste0("chr", chrom)]


### Methylation status of DMRs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data <- meth |>
  group_by(type) |>
  summarise(N = n()) |>
  mutate(
    percent = round(N / sum(N) * 100, 1),
    type = str_to_title(type),
    type = fct_inorder(type),
    csum = rev(cumsum(rev(N))),
    pos = N / 2 + lead(csum, 1), 
    pos = if_else(is.na(pos), N / 2, pos))

col <- brewer.pal(4, "Set2")[c(4, 3)]
meth_status <- ggplot(data, aes(x = "", y = N, fill = type)) +
  geom_col(width = 1, color = "black") +
  coord_polar(theta = "y", direction = -1) +
  labs(fill = "", title = "Methylation status (DMRs)") +
  theme_void() +
  theme(
    legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom") +
  scale_fill_manual(values = col) +
  geom_label_repel(aes(
    y = pos,
    label = paste0(comma(N), " (", percent, "%)")),
    size = 4,
    show.legend = F,
    nudge_x = .5,
    min.segment.length = 0,
    seed = 1337) +
  guides(fill = guide_legend(ncol = 1, keywidth = 1, keyheight = 1))


### Genomic regions of DMRs  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

col <- c(brewer.pal(4, "Set3"), "grey60")

data <- meth %>%
  group_by(region) %>%
  summarise(N = n()) %>%
  arrange(factor(region, levels = c("promoter", "exon", "intron", "intragenic", "intergenic"))) %>%
  mutate(percent = round(N / sum(N) * 100, 1),
         region = str_to_title(region),
         region = fct_inorder(region),
         csum = rev(cumsum(rev(N))),
         pos = N / 2 + lead(csum, 1), 
         pos = if_else(is.na(pos), N / 2, pos),
         col = col)

genomic_regions <- ggplot(data, aes(x = "", y = N, fill = region)) +
  geom_col(color = "black") +
  coord_polar(theta = "y", direction = -1) +
  labs(fill = "", title = "Genomic region (DMRs)") +
  theme_void() +
  theme(legend.text = element_text(size = 11), 
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = data$col) +
  geom_label_repel(aes(y = pos, label = paste0(comma(N), " (", percent, "%)")), 
                   size = 4,
                   show.legend = F,
                   nudge_x = .5,
                   min.segment.length = 0,
                   seed = 1337) +
  guides(fill = guide_legend(ncol = 2, nrow = 4, keywidth = 1, keyheight = 1))


### DM CpG-islands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Make a color gradient between 2 colors
color_fun <- function(col1, col2, n, plot = FALSE) {
  fc <- colorRampPalette(c(col1, col2))

  if (plot) pie(rep(1,n), col = fc(n))

  return(fc(n))
}

data <- data.frame(
    type = c("DM CpG-islands", "Not DM"),
    N = c(sum(!is.na(unique(meth$cgi))), length(cpg_islands))) %>%
  mutate(percent = round(N / sum(N) * 100, 1),
    type = fct_inorder(type),
    csum = rev(cumsum(rev(N))),
    pos = N / 2 + lead(csum, 1), 
    pos = if_else(is.na(pos), N / 2, pos),
    col = color_fun("#9999ff", "#ffcc66", 2))

CGIs <- ggplot(data, aes(x = "", y = N, fill = type)) +
  geom_col(color = "black") +
  coord_polar(theta = "y", direction = -1) +
  labs(fill = "", title = "CpG-islands") +
  theme_void() +
  theme(legend.text = element_text(size = 11), 
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = data$col) +
  geom_label_repel(aes(y = pos, label = paste0(comma(N), " (", percent, "%)")), 
                   size = 4,
                   show.legend = F,
                   nudge_x = .5,
                   min.segment.length = 0,
                   seed = 1337) +
  guides(fill = guide_legend(ncol = 1, keywidth = 1, keyheight = 1))


### Genes overlapping DMRs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gene_ids <- subset(meth, region %in% c("promoter", "exon") | !is.na(cgi_id) & !is.na(gene_id))$gene_id
gene_ids <- sapply(gene_ids, function(x) {
  strsplit(x, "[,]")[[1]]
}) |> unlist() |> unname() |> unique()

col <- brewer.pal(3, "Set2")[c(1,3)]

data <- data.frame(
  type = c("DMG", "Not DM"),
  N = c(length(gene_ids), nrow(genes) - length(gene_ids))) %>%
  mutate(percent = round(N / sum(N) * 100, 1),
         type = fct_inorder(type),
         csum = rev(cumsum(rev(N))),
         pos = N / 2 + lead(csum, 1), 
         pos = if_else(is.na(pos), N / 2, pos),
         col = col)

DMGs <- ggplot(data, aes(x = "", y = N, fill = type)) +
  geom_col(color = "black") +
  coord_polar(theta = "y", direction = -1) +
  labs(fill = "", title = "Annotated genes") +
  theme_void() +
  theme(legend.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = data$col) +
  geom_label_repel(aes(y = pos, label = paste0(comma(N), " (", percent, "%)")), 
                   size = 4,
                   show.legend = F,
                   nudge_x = .5,
                   min.segment.length = 0,
                   seed = 1337) +
  guides(fill = guide_legend(ncol = 1, keywidth = 1, keyheight = 1))


### Gene types of DMRs overlapping with promoter, exon and/or CGI ~~~~~~~~~~~~~

gene_type <- sapply(gene_ids, function(x) {
  subset(genes, gene_id %in% x)$gene_type2
})

# Make names look better
gene_type <- sub("protein_coding", "Protein coding", gene_type)
gene_type <- sub("pseudogene", "Pseudogene", gene_type)

data <- data.frame(table(gene_type)) %>%
  rename(Freq = "N", gene_type = "type") %>%
  arrange(desc(N)) %>%
  mutate(type = fct_inorder(type),
         percent = round(N / sum(N) * 100, 1),
         type = fct_inorder(type),
         csum = rev(cumsum(rev(N))),
         pos = N / 2 + lead(csum, 1), 
         pos = if_else(is.na(pos), N / 2, pos),
         col = brewer.pal(8,"Set2")[1:5])

gene_types <- ggplot(data, aes(x = "", y = N, fill = type)) +
  geom_col(color = "black") +
  coord_polar(theta = "y", direction = -1) +
  labs(fill = "", title = "Gene types (DMGs)") +
  theme_void() +
  theme(legend.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = data$col) +
  geom_label_repel(aes(y = pos, label = paste0(comma(N), " (", percent, "%)")), 
                   size = 4,
                   show.legend = F,
                   nudge_x = .5,
                   min.segment.length = 0,
                   seed = 1337) +
  guides(fill = guide_legend(ncol = 2, keywidth = 1, keyheight = 1))


### Save figure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

p <- plot_grid(
  meth_status,
  genomic_regions,
  CGIs,
  DMGs,
  gene_types,
  labels = "AUTO",
  align = "h"
)

filename <- "figure1.png"

ggsave(
  filename = filename,
  path = "images",
  plot = p,
  width = 9,
  height = 9,
  bg = "white",
  create.dir = TRUE
)