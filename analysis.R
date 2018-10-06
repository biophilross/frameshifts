#
# Script to perform Gene Ontology Enrichment analysis and GTEX tissue expression
# analysis.
#

# Necessary libraries in order to do this analysis

library(tidyverse)
library(cowplot)
library(org.Hs.eg.db)
library(clusterProfiler)

# Variables
hits_file <- "~/Downloads/conservative_reciprocal_best_hits.tsv"
gtex_file <- "~/Desktop/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct"
tpm_threshold <- 2

# GO enrichment

hits <- readr::read_tsv(hits_file, col_names = T)

genes <- hits$query

# we need the EntrzIDs of these genes to use enrichGO
# use the biological ID translator function for this
gene.df <- clusterProfiler::bitr(genes, fromType = "ENSEMBL",
                toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

gobp <- clusterProfiler::enrichGO(
                gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.10)

gomf <- clusterProfiler::enrichGO(
                 gene          = gene.df$ENTREZID,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.10)

# Plot these enrichments
g <- clusterProfiler::dotplot(gobp, showCategory=30)
cowplot::save_plot("biological_process_go_enrichment.png", plot=g)

g <- clusterProfiler::dotplot(gomf, showCategory=30)
cowplot::save_plot("molecular_function_go_enrichment.png", plot=g)

# GTEX tissue expression analysis

gtex <- readr::read_tsv(gtex_file, skip = 2, col_names = T) %>%
  dplyr::select(-Description) %>%
  tidyr::gather(key = tissue, value = exp, -gene_id)

# remove annoying ".#" parts of each gene ID
gtex$gene_id <- sapply(stringr::str_split(gtex$gene_id,"[.]"), function(x) {x[1]})

# only include genes in our gene list
fgtex <- gtex %>% filter(gene_id %in% genes)

# group by gene ID to create histogram of how many tissues these genes are
# expressed in
gene_expression <- fgtex %>% 
  dplyr::group_by(gene_id) %>% 
  dplyr::summarise(n = sum(exp > tpm_threshold))

g <- ggplot(gene_expression, aes(x=n)) + 
  geom_histogram(color="grey70") +
  xlab("Number of Tissues") + 
  ylab("Number of Genes")

cowplot::save_plot("gene_expression_distribution.png", plot=g)

# group by tissue to create a table of how many genes are express in each 
# tissue type
tissue_expression <- fgtex %>% 
  dplyr::filter(exp > tpm_threshold) %>% 
  dplyr::group_by(tissue) %>% 
  dplyr::summarise(n = n()) %>%
  dplyr::arrange(desc(n))

g <- ggplot(tissue_expression, aes(x=reorder(tissue, n),y=n)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  ylab("Number of genes") + 
  xlab("")

cowplot::save_plot("tissue_expression_barplot.png", plot=g)