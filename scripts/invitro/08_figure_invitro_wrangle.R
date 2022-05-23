# Author: Emilano Sotelo & Melyssa Minto


# Load packages -----------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)

# Read in data ------------------------------------------------------------
# Read differentially expressed genes 
knockdown_genes_zic1 <- read_tsv("../../results/invitro/diffexpr/res_zic1_ashr.tsv")
knockdown_genes_zic2 <- read_tsv("../../results/invitro/diffexpr/res_zic2_ashr.tsv")
dev_genes <- read_tsv("../../results/invitro/diffexpr/res_d7_d3_ashr.tsv")

# Read In vitro Zic data
zic_invitro <- read_tsv("../../results/FinalTables/invitro_mapped_data.txt")
zic_invivo <- read_tsv("../../results/FinalTables/mapped_data.txt")

# Zic diff info
zic_diff_peaks <- as.data.frame(readRDS("../../results/invitro/diffbind_cutnrun_zic/peaks_7v3_all_annotation.Rds"))


# Read RRHO gene sets
zic1_up_rrho_genes <- read_csv("../../results/invitro/rrho/RRHO_GO_MostUpregulatedZic1_KD_Ctrl_VS_D7_D3.csv",col_names = FALSE)
zic1_down_rrho_genes <- read_csv("../../results/invitro/rrho/RRHO_GO_MostDownregulatedZic1_KD_Ctrl_VS_D7_D3.csv",col_names = FALSE)

zic2_up_rrho_genes <- read_csv("../../results/invitro/rrho/RRHO_GO_MostUpregulatedZic2_KD_Ctrl_VS_D7_D3.csv",col_names = FALSE)
zic2_down_rrho_genes <- read_csv("../../results/invitro/rrho/RRHO_GO_MostDownregulatedZic2_KD_Ctrl_VS_D7_D3.csv",col_names = FALSE)



# Data Wrangling ----------------------------------------------------------
# adding regulation and peak identifier
zic_diff_peaks <- zic_diff_peaks %>%
  dplyr::mutate(zic_sig = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "UP",
    padj < 0.05 & log2FoldChange < 0 ~ "DOWN",
    padj > 0.05 ~ "N.S.",
    is.na(padj) ~ "NA"
  ), zic_peak = paste0(seqnames,":",start - 1,"-",end)) %>%
  dplyr::select(zic_peak, zic_sig)

# how many peaks overalpped loops
summary(zic_diff_peaks$zic_peak %in% zic_invitro$zic_peak)

# filtering for peaks within loops
zic_diff_peaks <- zic_diff_peaks %>% 
  dplyr::filter(zic_peak %in% zic_invitro$zic_peak)



# > get number of peaks for gene ------------------------------------------
# Get number of peaks in vivo and in vitro and by direction
peaks_per_gene <- zic_invivo %>%
  dplyr::select(gene_name, zic_peak) %>%
  dplyr::distinct() %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(npeaks = n()) %>% 
  ungroup() 

peaks_per_gene_dir <- zic_invivo %>%
  dplyr::select(gene_name, zic_peak,zic_sig) %>%
  dplyr::distinct() %>%
  dplyr::group_by(gene_name, zic_sig) %>%
  dplyr::summarise(npeaks = n()) %>% 
  ungroup() 

peaks_per_gene_invitro <- zic_invitro %>%
  dplyr::select(gene_name, zic_peak) %>%
  dplyr::distinct() %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(npeaks = n()) %>% 
  ungroup() 

peaks_per_gene_invitro_dir <- zic_invitro %>%
  left_join(zic_diff_peaks, by = "zic_peak") %>%
  dplyr::select(gene_name, zic_peak,zic_sig) %>%
  dplyr::distinct() %>%
  dplyr::group_by(gene_name, zic_sig) %>%
  dplyr::summarise(npeaks = n()) %>% 
  ungroup()


colnames(peaks_per_gene_invitro_dir) <- paste0(colnames(peaks_per_gene_invitro_dir), "_invitro")


#> merge differential expression data --------------------------------------

# Change column names before the join 
colnames(knockdown_genes_zic1) <- paste0("zic1_",colnames(knockdown_genes_zic1))
colnames(knockdown_genes_zic2) <- paste0("zic2_",colnames(knockdown_genes_zic2))
colnames(dev_genes) <- paste0("dev_",colnames(dev_genes))

# join peaks and gene expression data
knockdown_genes_zic1_npeaks <- knockdown_genes_zic1 %>% 
  # adding peak-gene mapping
  left_join( peaks_per_gene, by = c("zic1_gene_names" = "gene_name" )) %>% 
  # addong gene expresison
  left_join(dev_genes, by = c("zic1_gene_names" = "dev_gene_names")) %>%
  dplyr::mutate(opposite_direction = zic1_log2FoldChange * dev_log2FoldChange < 0) 

knockdown_genes_zic2_npeaks <- knockdown_genes_zic2 %>% 
  left_join( peaks_per_gene, by = c("zic2_gene_names" = "gene_name" )) %>% 
  left_join( dev_genes, by = c("zic2_gene_names" = "dev_gene_names")) %>%
  dplyr::mutate(opposite_direction = zic2_log2FoldChange * dev_log2FoldChange < 0)

knockdown_genes_zic1_npeaks_in_vitro <- knockdown_genes_zic1 %>% 
  left_join(peaks_per_gene_invitro, by = c("zic1_gene_names" = "gene_name" )) %>% 
  left_join(dev_genes, by = c("zic1_gene_names" = "dev_gene_names")) %>%
  dplyr::mutate(opposite_direction = zic1_log2FoldChange * dev_log2FoldChange < 0)

knockdown_genes_zic2_npeaks_in_vitro <- knockdown_genes_zic2 %>% 
  left_join(peaks_per_gene_invitro, by = c("zic2_gene_names" = "gene_name" )) %>% 
  left_join(dev_genes, by = c("zic2_gene_names" = "dev_gene_names")) %>%
  dplyr::mutate(opposite_direction = zic2_log2FoldChange * dev_log2FoldChange < 0)

knockdown_genes_zic1_npeaks_dir <- knockdown_genes_zic1 %>% 
  left_join(peaks_per_gene_dir, by = c("zic1_gene_names" = "gene_name" )) %>% 
  left_join(dev_genes, by = c("zic1_gene_names" = "dev_gene_names")) %>%
  dplyr::mutate(opposite_direction = zic1_log2FoldChange * dev_log2FoldChange < 0)

knockdown_genes_zic2_npeaks_dir <-knockdown_genes_zic2 %>% 
  left_join(peaks_per_gene_dir, by = c("zic2_gene_names" = "gene_name" )) %>% 
  left_join(dev_genes, by = c("zic2_gene_names" = "dev_gene_names")) %>%
  dplyr::mutate(opposite_direction = zic2_log2FoldChange * dev_log2FoldChange < 0)

knockdown_genes_zic1_npeaks_in_vitro_dir  <- knockdown_genes_zic1 %>% 
  left_join(peaks_per_gene_invitro_dir, by = c("zic1_gene_names" = "gene_name_invitro" )) %>% 
  left_join(dev_genes, by = c("zic1_gene_names" = "dev_gene_names")) %>%
  dplyr::mutate(opposite_direction = zic1_log2FoldChange * dev_log2FoldChange < 0)

knockdown_genes_zic2_npeaks_in_vitro_dir <-  knockdown_genes_zic2 %>%
  left_join(peaks_per_gene_invitro_dir, by = c("zic2_gene_names" = "gene_name_invitro" )) %>% 
  left_join(dev_genes, by = c("zic2_gene_names" = "dev_gene_names")) %>%
  dplyr::mutate(opposite_direction = zic2_log2FoldChange * dev_log2FoldChange < 0)


# > Add RRHO data ---------------------------------------------------------
knockdown_genes_zic1_npeaks <- knockdown_genes_zic1_npeaks %>%
  dplyr::mutate(rrho_category = case_when(
    zic1_gene_names %in% zic1_up_rrho_genes$X1 ~ "up",
    zic1_gene_names %in% zic1_down_rrho_genes$X1 ~ "down",
    TRUE ~ "ns"
  ))

knockdown_genes_zic2_npeaks <- knockdown_genes_zic2_npeaks %>%
  dplyr::mutate(rrho_category = case_when(
    zic2_gene_names %in% zic2_up_rrho_genes$X1 ~ "up",
    zic2_gene_names %in% zic2_down_rrho_genes$X1 ~ "down",
    TRUE ~ "ns"
  ))


# > Get Zic dependent and independent genes -------------------------------
peaks_per_gene_invitro_dir <- peaks_per_gene_invitro_dir %>%
  dplyr::filter(zic_sig_invitro != "NA" & !is.na(zic_sig_invitro))

knockdown_genes_zic1_npeaks_dir_zic_dependent <- knockdown_genes_zic1_npeaks_dir %>%
  left_join(peaks_per_gene_invitro_dir, by  = c("zic1_gene_names" = "gene_name_invitro",
                                                "zic_sig"="zic_sig_invitro")) %>%
  drop_na(contains("padj")) %>%
  dplyr::filter(zic1_padj < 0.05 & dev_padj < 0.05) %>%
  dplyr::mutate(diff_category = case_when(
    zic1_log2FoldChange < 0 & dev_log2FoldChange > 0  ~ "dependent_d7",
    zic1_log2FoldChange > 0 & dev_log2FoldChange < 0  ~ "dependent_d3",
    zic1_log2FoldChange > 0 & dev_log2FoldChange > 0  ~ "independent_d7",
    zic1_log2FoldChange < 0 & dev_log2FoldChange < 0  ~ "independent_d3",
  ))

knockdown_genes_zic2_npeaks_dir_zic_dependent <- knockdown_genes_zic2_npeaks_dir %>%
  left_join(peaks_per_gene_invitro_dir, by  = c("zic2_gene_names" = "gene_name_invitro",
                                                "zic_sig"="zic_sig_invitro")) %>%
  drop_na(contains("padj")) %>%
  dplyr::filter(zic2_padj < 0.05 & dev_padj < 0.05) %>%
  dplyr::mutate(diff_category = case_when(
    zic2_log2FoldChange < 0 & dev_log2FoldChange > 0  ~ "dependent_d7",
    zic2_log2FoldChange > 0 & dev_log2FoldChange < 0  ~ "dependent_d3",
    zic2_log2FoldChange > 0 & dev_log2FoldChange > 0  ~ "independent_d7",
    zic2_log2FoldChange < 0 & dev_log2FoldChange < 0  ~ "independent_d3",
  )) 



# > Adding GO enrichment --------------------------------------------------

# Wrangling for go enrichment
zic1_genes_by_category <- knockdown_genes_zic1_npeaks_dir_zic_dependent %>%
  dplyr::select(zic1_gene_names, diff_category)

zic1_genes_by_category <- split(zic1_genes_by_category$zic1_gene_names, zic1_genes_by_category$diff_category)

zic2_genes_by_category <- knockdown_genes_zic2_npeaks_dir_zic_dependent %>%
  dplyr::select(zic2_gene_names, diff_category)

zic2_genes_by_category <- split(zic2_genes_by_category$zic2_gene_names, zic2_genes_by_category$diff_category)

##### GO enrichments

zic1_comp_go <- compareCluster(zic1_genes_by_category,
                               OrgDb = org.Mm.eg.db,
                               fun           = "enrichGO",
                               keyType = "SYMBOL",
                               pvalueCutoff  = 0.05,
                               pAdjustMethod = "BH",
                               ont = "BP")

zic2_comp_go <- compareCluster(zic2_genes_by_category,
                               OrgDb = org.Mm.eg.db,
                               fun           = "enrichGO",
                               keyType = "SYMBOL",
                               pvalueCutoff  = 0.05,
                               pAdjustMethod = "BH",
                               ont = "BP")



# Output results ----------------------------------------------------------
# Npeaks data
write_tsv(knockdown_genes_zic1_npeaks,"../../results/invitro/figure_invitro/knockdown_genes_zic1_npeaks.tsv")
write_tsv(knockdown_genes_zic2_npeaks,"../../results/invitro/figure_invitro/knockdown_genes_zic2_npeaks.tsv")

write_tsv(knockdown_genes_zic1_npeaks_in_vitro,"../../results/invitro/figure_invitro/knockdown_genes_zic1_npeaks_in_vitro.tsv")
write_tsv(knockdown_genes_zic2_npeaks_in_vitro,"../../results/invitro/figure_invitro/knockdown_genes_zic2_npeaks_in_vitro")

write_tsv(knockdown_genes_zic1_npeaks_dir,"../../results/invitro/figure_invitro/knockdown_genes_zic1_npeaks_dir.tsv")
write_tsv(knockdown_genes_zic2_npeaks_dir,"../../results/invitro/figure_invitro/knockdown_genes_zic2_npeaks_dir.tsv")

write_tsv(knockdown_genes_zic1_npeaks_in_vitro_dir,"../../results/invitro/figure_invitro/knockdown_genes_zic1_npeaks_in_vitro_dir.tsv")
write_tsv(knockdown_genes_zic2_npeaks_in_vitro_dir,"../../results/invitro/figure_invitro/knockdown_genes_zic2_npeaks_in_vitro_dir.tsv")

write_tsv(knockdown_genes_zic1_npeaks_dir_zic_dependent,"../../results/invitro/figure_invitro/knockdown_genes_zic1_npeaks_in_vitro_dir.tsv")
write_tsv(knockdown_genes_zic2_npeaks_dir_zic_dependent,"../../results/invitro/figure_invitro/knockdown_genes_zic2_npeaks_in_vitro_dir.tsv")


# GO data
saveRDS(zic1_comp_go, "../../results/invitro/figure_invitro/zic1_comp_go.Rds")
saveRDS(zic2_comp_go, "../../results/invitro/figure_invitro/zic2_comp_go.Rds")

