# Author: Emiliano Sotelo & Melyssa Minto
# This script will perform a rank-rank hypergeometric overlap test on Zic1 and Zic2 KD RNA-seq results


# Load libraries ----------------------------------------------------------
library(tidyverse)
library(RRHO)


# Read in data ------------------------------------------------------------
# Read differentially expressed genes 
knockdown_genes_zic1 <- read_tsv("../../results/invitro/diffexpr/res_zic1_ashr.tsv")
knockdown_genes_zic2 <- read_tsv("../../results/invitro/diffexpr/res_zic2_ashr.tsv")
dev_genes <- read_tsv("../../results/invitro/diffexpr/res_d7_d3_ashr.tsv")



# Wrangle data ------------------------------------------------------------


# > Merging differential analysis results ---------------------------------

# Change column names before the join 
colnames(knockdown_genes_zic1) <- paste0("zic1_",colnames(knockdown_genes_zic1))
colnames(knockdown_genes_zic2) <- paste0("zic2_",colnames(knockdown_genes_zic2))
colnames(dev_genes) <- paste0("dev_",colnames(dev_genes))

# Join results tables
knockdown_genes_zic1_dev <- inner_join(knockdown_genes_zic1,
                                       dev_genes, 
                                        by = c("zic1_gene_names" = "dev_gene_names"), 
                                        na_matches = "never")
knockdown_genes_zic1_dev <- knockdown_genes_zic1_dev[!duplicated(knockdown_genes_zic1_dev$zic1_gene_names) | duplicated(knockdown_genes_zic1_dev$zic1_gene_names,fromLast = TRUE),]

knockdown_genes_zic2_dev <- inner_join(knockdown_genes_zic2,
                                        dev_genes, 
                                        by = c("zic2_gene_names" = "dev_gene_names"),
                                        na_matches = "never")
knockdown_genes_zic2_dev <- knockdown_genes_zic2_dev[!duplicated(knockdown_genes_zic2_dev$zic2_gene_names) | duplicated(knockdown_genes_zic2_dev$zic2_gene_names,fromLast = TRUE),]


# > Prepare tables for rrho -----------------------------------------------
signal_zic1_vs_dev <- knockdown_genes_zic1_dev %>% 
  select(contains("adj"),contains("log2"),zic1_gene_names) %>%
  distinct() %>%
  drop_na() %>%
  mutate(sgnlogpadj = log10(zic1_padj)* sign(as.numeric(zic1_log2FoldChange))) %>% 
  dplyr::rename( SYMBOL = zic1_gene_names) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

signal_dev_vs_zic1_neg <- knockdown_genes_zic1_dev %>% 
  select(contains("adj"),contains("log2"),zic1_gene_names) %>%
  distinct() %>%
  drop_na() %>%
  mutate(sgnlogpadj = - log10(dev_padj)* sign(as.numeric(dev_log2FoldChange))) %>% 
  dplyr::rename( SYMBOL = zic1_gene_names) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

signal_zic2_vs_dev <- knockdown_genes_zic2_dev %>% 
  select(contains("adj"),contains("log2"),zic2_gene_names) %>%
  distinct() %>%
  drop_na() %>%
  mutate(sgnlogpadj = log10(zic2_padj)* sign(as.numeric(zic2_log2FoldChange))) %>% 
  dplyr::rename( SYMBOL = zic2_gene_names) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

signal_dev_vs_zic2_neg <- knockdown_genes_zic2_dev %>% 
  select(contains("adj"),contains("log2"),zic2_gene_names) %>%
  distinct() %>%
  drop_na() %>%
  mutate(sgnlogpadj = - log10(dev_padj)* sign(as.numeric(dev_log2FoldChange))) %>% 
  dplyr::rename( SYMBOL = zic2_gene_names) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()


# Run RRHO analysis -------------------------------------------------------
RRHO(as.data.frame(signal_zic1_vs_dev),
                as.data.frame(signal_dev_vs_zic1_neg),
                BY=TRUE,
                stepsize = 100,
                alternative='two.sided',
                log10.ind = TRUE,
                plot = TRUE,
                outputdir = "../../results/invitro/rrho/",
                labels = c("Zic1_KD_Ctrl","D7_D3"))

RRHO(as.data.frame(signal_zic2_vs_dev),
                as.data.frame(signal_dev_vs_zic2_neg),
                BY=TRUE,
                stepsize = 100,
                alternative='two.sided',
                log10.ind = TRUE,
                plot = TRUE,
                outputdir = "../../results/invitro/rrho/",
                labels = c("Zic2_KD_Ctrl","D7_D3"))



