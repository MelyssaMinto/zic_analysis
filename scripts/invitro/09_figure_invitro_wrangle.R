# Author: Emilano Sotelo & Melyssa Minto


# Load packages -----------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ChIPseeker)
# Read in data ------------------------------------------------------------
# Read differentially expressed genes 
knockdown_genes_zic1 <- read_tsv("../../results/invitro/diffexpr/res_zic1_ashr.tsv")
knockdown_genes_zic2 <- read_tsv("../../results/invitro/diffexpr/res_zic2_ashr.tsv")
dev_genes <- read_tsv("../../results/invitro/diffexpr/res_d7_d3_ashr.tsv")

# Read In vitro Zic data
zic_invitro <- read_tsv("../../results/FinalTables/invitro_mapped_data.txt")
zic_invivo <- read_tsv("../../results/FinalTables/mapped_data.txt")
loop_data <- read_tsv("../../results/FinalTables/loop_data.txt")

# Zic diff info
zic_diff_peaks <- as.data.frame(readRDS("../../results/invitro/diffbind_cutnrun_zic/peaks_7v3_all_annotation.Rds"))


# Read RRHO gene sets
# > Genes that are concordantly upregulated in Zic1 KD and Enriched in DIV3
zic1_up_rrho_genes <- read_csv("../../results/invitro/rrho/RRHO_GO_MostUpregulatedZic1KD v WT_VS_D3_D7.csv",col_names = FALSE)
# > Genes that are concordantly downregulated in Zic1 KD and Enriched in DIV7
zic1_down_rrho_genes <- read_csv("../../results/invitro/rrho/RRHO_GO_MostDownregulatedZic1KD v WT_VS_D3_D7.csv",col_names = FALSE)
# >  Genes that are concordantly upregulated in Zic2 KD and Enriched in DIV3
zic2_up_rrho_genes <- read_csv("../../results/invitro/rrho/RRHO_GO_MostUpregulatedZic2KD v WT_VS_D3_D7.csv",col_names = FALSE)
#> Genes that are concordantly upregulated in Zic2 KD and Enriched in DIV7
zic2_down_rrho_genes <- read_csv("../../results/invitro/rrho/RRHO_GO_MostDownregulatedZic2KD v WT_VS_D3_D7.csv",col_names = FALSE)


# Data Wrangling ----------------------------------------------------------
# adding regulation and peak identifier
zic_diff_peaks <- zic_diff_peaks %>%
  dplyr::mutate(zic_sig = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "UP",
    padj < 0.05 & log2FoldChange < 0 ~ "DOWN",
    padj > 0.05 ~ "N.S.",
    is.na(padj) ~ "NA"
  ), zic_peak = paste0(seqnames,":",start - 1,"-",end)) %>%
  dplyr::select(zic_peak, zic_peak_sig = zic_sig, zic_peak_baseMean =  baseMean, zic_peak_log2FoldChange = log2FoldChange, zic_peak_lfcSE = lfcSE, zic_peak_pvalue = pvalue, zic_peak_padj = padj )

# how many peaks overalpped loops
summary(zic_diff_peaks$zic_peak %in% zic_invitro$zic_peak)

# filtering for peaks within loops
# zic_diff_peaks <- zic_diff_peaks %>% 
#   dplyr::filter(zic_peak %in% zic_invitro$zic_peak)



#> merge differential expression data --------------------------------------

# Change column names before the join 
colnames(knockdown_genes_zic1) <- paste0("zic1_",colnames(knockdown_genes_zic1))
colnames(knockdown_genes_zic2) <- paste0("zic2_",colnames(knockdown_genes_zic2))
colnames(dev_genes) <- paste0("dev_",colnames(dev_genes))

# join peaks and gene expression data


knockdown_genes_zic1_npeaks_in_vitro <- knockdown_genes_zic1 %>% 
  left_join(dev_genes, by = c("zic1_gene_names" = "dev_gene_names")) %>%
  dplyr::mutate(opposite_direction = zic1_log2FoldChange * dev_log2FoldChange < 0)

knockdown_genes_zic2_npeaks_in_vitro <- knockdown_genes_zic2 %>% 
  left_join(dev_genes, by = c("zic2_gene_names" = "dev_gene_names")) %>%
  dplyr::mutate(opposite_direction = zic2_log2FoldChange * dev_log2FoldChange < 0)


# > Add RRHO data ---------------------------------------------------------
knockdown_genes_zic1_npeaks_in_vitro <- knockdown_genes_zic1_npeaks_in_vitro %>%
  dplyr::mutate(rrho_category = case_when(
    zic1_gene_names %in% zic1_up_rrho_genes$X1 ~ "up", # early genes that remain on
    zic1_gene_names %in% zic1_down_rrho_genes$X1 ~ "down", #late genes that don't upregulate
    TRUE ~ "ns"
  ))

knockdown_genes_zic2_npeaks_in_vitro <- knockdown_genes_zic2_npeaks_in_vitro %>%
  dplyr::mutate(rrho_category = case_when(
    zic2_gene_names %in% zic2_up_rrho_genes$X1 ~ "up", # early genes that remain on
    zic2_gene_names %in% zic2_down_rrho_genes$X1 ~ "down", # late genes that don't upregulate
    TRUE ~ "ns"
  ))


# > Get Zic dependent and independent genes -------------------------------
# peaks_per_gene_invitro_dir <- peaks_per_gene_invitro_dir %>%
#   dplyr::filter(zic_sig_invitro != "NA" & !is.na(zic_sig_invitro))
# 
# knockdown_genes_zic1_npeaks_dir_zic_dependent <- knockdown_genes_zic1_npeaks_dir %>%
#   left_join(peaks_per_gene_invitro_dir, by  = c("zic1_gene_names" = "gene_name_invitro",
#                                                 "zic_sig"="zic_sig_invitro")) %>%
#   drop_na(contains("padj")) %>%
#   dplyr::filter(zic1_padj < 0.05 & dev_padj < 0.05) %>%
#   dplyr::mutate(diff_category = case_when(
#     zic1_log2FoldChange < 0 & dev_log2FoldChange > 0  ~ "dependent_d7",
#     zic1_log2FoldChange > 0 & dev_log2FoldChange < 0  ~ "dependent_d3",
#     zic1_log2FoldChange > 0 & dev_log2FoldChange > 0  ~ "independent_d7",
#     zic1_log2FoldChange < 0 & dev_log2FoldChange < 0  ~ "independent_d3",
#   ))
# 
# knockdown_genes_zic2_npeaks_dir_zic_dependent <- knockdown_genes_zic2_npeaks_dir %>%
#   left_join(peaks_per_gene_invitro_dir, by  = c("zic2_gene_names" = "gene_name_invitro",
#                                                 "zic_sig"="zic_sig_invitro")) %>%
#   drop_na(contains("padj")) %>%
#   dplyr::filter(zic2_padj < 0.05 & dev_padj < 0.05) %>%
#   dplyr::mutate(diff_category = case_when(
#     zic2_log2FoldChange < 0 & dev_log2FoldChange > 0  ~ "dependent_d7",
#     zic2_log2FoldChange > 0 & dev_log2FoldChange < 0  ~ "dependent_d3",
#     zic2_log2FoldChange > 0 & dev_log2FoldChange > 0  ~ "independent_d7",
#     zic2_log2FoldChange < 0 & dev_log2FoldChange < 0  ~ "independent_d3",
#   )) 



# > make invitro data -----------------------------------------------------

invitro_data =
  # combine KD data
  bind_rows(list( Zic1 = knockdown_genes_zic1_npeaks_in_vitro %>% setNames(tolower(gsub("zic1_","",names(.)))),
                               Zic2  = knockdown_genes_zic2_npeaks_in_vitro %>% setNames(tolower(gsub("zic2_","",names(.))))),
                         .id = "id") %>% 
  drop_na(contains("padj")) %>%
  dplyr::filter(rrho_category != "ns") %>% 
  # adding categories for genes
  dplyr::mutate(cat = case_when( padj < 0.05 & dev_padj < 0.05 ~ "Zic Dependent Developmental",
                                 padj >= 0.05 & dev_padj < 0.05 ~ "Developmental",
                                 padj < 0.05 & dev_padj >= 0.05 ~ "Zic Dependent",
                                 TRUE ~ "N.S."),
                cat = fct_relevel(cat, c("Developmental", "Zic Dependent", "Zic Dependent Developmental")),
                id = paste0(id, " KD")) %>% 
  dplyr::filter(cat != "N.S.") %>% 
  # determining whether a gene was affected by Zic1 KD, Zic2 KD, or both
  group_by(gene_names) %>% 
  dplyr::mutate(bothKD = case_when(n_distinct(id) > 1 ~ "Zic1 & Zic2 KD",
                                   n_distinct(id) == 1 ~ id)  ) %>% 
  ungroup() %>% 
  dplyr::mutate(label = ifelse(bothKD == "Zic1 & Zic2", gene_names, "")) %>% 
  # adding gene_peak mapping
  left_join(zic_invitro %>% dplyr::rename( loop_id = id), by = c("gene_names" = "gene_name")) %>% 
  # adding differential zic binding data
  left_join(zic_diff_peaks, by = "zic_peak" ) %>%  
  # count peaks and order by npeaks
  group_by(gene_names) %>% 
  dplyr::mutate(num_peaks = n_distinct(zic_peak)) %>%
  dplyr::mutate(num_peaks = ifelse( n() == sum(is.na(zic_peak)), 0, num_peaks) ) %>% 
  ungroup() %>% 
  # create differential cut offs
  dplyr::mutate(dev_sig = case_when( dev_log2foldchange > 1 & dev_padj < 0.05 ~ "Up",
                                     dev_log2foldchange < -1 & dev_padj < 0.05 ~ "Down",
                                     TRUE ~ "N.S.")) %>% 
  dplyr::mutate(kd_sig = case_when( log2foldchange > 1 & padj < 0.05 ~ "Up",
                                    log2foldchange < -1 & padj < 0.05 ~ "Down",
                                    TRUE ~ "N.S.")) %>% 
  dplyr::mutate(zic_sig = case_when( zic_peak_log2FoldChange > 0 & zic_peak_padj < 0.05 ~ "Up",
                                     zic_peak_log2FoldChange < 0 & zic_peak_padj < 0.05 ~ "Down",
                                     TRUE ~ "Static")) 
 

# > export Dev, Zic Dep, and ZDD gene regulatory regions ------------------------------------

invitro_data %>% 
  dplyr::filter(cat == "Zic Dependent Developmental") %>% 
  dplyr::select(loop_id) %>% 
  dplyr::left_join(loop_data, by = c("loop_id" = "id")) %>% 
  dplyr::select(starts_with("anchor")) %>% 
  pivot_longer(cols = starts_with("anchor")) %>% 
  dplyr::mutate(chr = str_extract(value,"[^\\:]+") ,
                start = str_extract(value, "(?<=:).+(?=-)") , 
                end = str_extract(value, "\\b\\w+$") ) %>% 
  dplyr::select(-name, -value) %>% 
  distinct() %>% 
  drop_na() %>% 
  dplyr::mutate(name = ".", score = 1, strand = ".") %>% 
  write_tsv("../../results/invitro/beds/zdd_anchors.bed", col_names = F)

invitro_data %>% 
  dplyr::filter(cat == "Zic Dependent") %>% 
  dplyr::select(loop_id) %>% 
  dplyr::left_join(loop_data, by = c("loop_id" = "id")) %>% 
  dplyr::select(starts_with("anchor")) %>% 
  pivot_longer(cols = starts_with("anchor")) %>% 
  dplyr::mutate(chr = str_extract(value,"[^\\:]+") ,
                start = str_extract(value, "(?<=:).+(?=-)") , 
                end = str_extract(value, "\\b\\w+$") ) %>% 
  dplyr::select(-name, -value) %>% 
  distinct() %>% 
  drop_na() %>% 
  dplyr::mutate(name = ".", score = 1, strand = ".") %>% 
  write_tsv("../../results/invitro/beds/zd_anchors.bed", col_names = F)

invitro_data %>% 
  dplyr::filter(cat == "Developmental") %>% 
  dplyr::select(loop_id) %>% 
  dplyr::left_join(loop_data, by = c("loop_id" = "id")) %>% 
  dplyr::select(starts_with("anchor")) %>% 
  pivot_longer(cols = starts_with("anchor")) %>% 
  dplyr::mutate(chr = str_extract(value,"[^\\:]+") ,
                start = str_extract(value, "(?<=:).+(?=-)") , 
                end = str_extract(value, "\\b\\w+$") ) %>% 
  dplyr::select(-name, -value) %>% 
  distinct() %>% 
  drop_na() %>% 
  dplyr::mutate(name = ".", score = 1, strand = ".") %>% 
  write_tsv("../../results/invitro/beds/dev_anchors.bed", col_names = F)


invitro_data %>% 
  dplyr::filter(cat == "Zic Dependent Developmental") %>% 
  dplyr::select(zic_peak) %>% 
  dplyr::mutate(chr = str_extract(zic_peak,"[^\\:]+") ,
                start = str_extract(zic_peak, "(?<=:).+(?=-)") , 
                end = str_extract(zic_peak, "\\b\\w+$") ) %>% 
  dplyr::select(-zic_peak) %>% 
  distinct() %>% 
  drop_na() %>% 
  dplyr::mutate(name = ".", score = 1, strand = ".") %>% 
  write_tsv("../../results/invitro/beds/zdd_zicPeaks.bed", col_names = F)

invitro_data %>% 
  dplyr::filter(cat == "Zic Dependent") %>% 
  dplyr::select(zic_peak) %>% 
  dplyr::mutate(chr = str_extract(zic_peak,"[^\\:]+") ,
                start = str_extract(zic_peak, "(?<=:).+(?=-)") , 
                end = str_extract(zic_peak, "\\b\\w+$") ) %>% 
  dplyr::select(-zic_peak) %>% 
  distinct() %>% 
  drop_na() %>% 
  dplyr::mutate(name = ".", score = 1, strand = ".") %>% 
  write_tsv("../../results/invitro/beds/zd_zicPeaks.bed", col_names = F)

invitro_data %>% 
  dplyr::filter(cat == "Developmental") %>% 
  dplyr::select(zic_peak) %>% 
  dplyr::mutate(chr = str_extract(zic_peak,"[^\\:]+") ,
                start = str_extract(zic_peak, "(?<=:).+(?=-)") , 
                end = str_extract(zic_peak, "\\b\\w+$") ) %>% 
  dplyr::select(-zic_peak) %>% 
  distinct() %>% 
  drop_na() %>% 
  dplyr::mutate(name = ".", score = 1, strand = ".") %>% 
  write_tsv("../../results/invitro/beds/dev_zicPeaks.bed", col_names = F)
  


# > Adding GO enrichment --------------------------------------------------
genes =
  invitro_data %>% 
  dplyr::select(gene_names) %>% 
  distinct() %>% 
  pull(gene_names)

gene_cnvr = clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

genes_cat = invitro_data %>% 
  left_join(gene_cnvr, by = c("gene_names" = "SYMBOL")) %>% 
  dplyr::filter(dev_sig != "N.S.") %>% 
  dplyr::select(gene_names,ENTREZID, cat, dev_sig, num_peaks) %>% 
  distinct()

genes_cat_targets = invitro_data %>% 
  left_join(gene_cnvr, by = c("gene_names" = "SYMBOL")) %>% 
  dplyr::filter(num_peaks > 0) %>% 
  dplyr::select(gene_names,ENTREZID, cat, dev_sig, num_peaks) %>% 
  distinct()

zdd_genes = genes_cat %>% 
  dplyr::filter(cat == "Zic Dependent Developmental") %>% 
  distinct() 

zdd_genes_targets = genes_cat %>% 
  dplyr::filter(cat == "Zic Dependent Developmental", num_peaks > 0) %>% 
  distinct() 

zdd_genes <- split(zdd_genes$gene_names, zdd_genes$dev_sig)
zdd_genes_targets <- split(zdd_genes_targets$gene_names, zdd_genes_targets$dev_sig)

g <- split(genes_cat_targets$gene_names, genes_cat_targets$cat)

g$N.S. <- NULL


# running GO enrichment 
g_bp <- compareCluster(g,
                       OrgDb = org.Mm.eg.db,
                       fun           = "enrichGO",
                       keyType = "SYMBOL",
                       pvalueCutoff  = 0.05,
                       pAdjustMethod = "BH",
                       ont = "BP")
g_mf <- compareCluster(g,
                       OrgDb = org.Mm.eg.db,
                       fun           = "enrichGO",
                       keyType = "SYMBOL",
                       pvalueCutoff  = 0.05,
                       pAdjustMethod = "BH",
                       ont = "MF")
g_cc <- compareCluster(g,
                       OrgDb = org.Mm.eg.db,
                       fun           = "enrichGO",
                       keyType = "SYMBOL",
                       pvalueCutoff  = 0.05,
                       pAdjustMethod = "BH",
                       ont = "CC")


g_cat_bp = compareCluster(gene_names ~ cat + dev_sig,
               data = genes_cat_targets,
               OrgDb = org.Mm.eg.db,
               fun           = "enrichGO",
               keyType = "SYMBOL",
               pvalueCutoff  = 0.05,
               pAdjustMethod = "BH",
               ont = "BP")

g_cat_mf = compareCluster(gene_names ~ cat + dev_sig,
                          data = genes_cat_targets,
                          OrgDb = org.Mm.eg.db,
                          fun           = "enrichGO",
                          keyType = "SYMBOL",
                          pvalueCutoff  = 0.05,
                          pAdjustMethod = "BH",
                          ont = "MF")



g_cat_cc = compareCluster(gene_names ~ cat + dev_sig,
                          data = genes_cat_targets,
                          OrgDb = org.Mm.eg.db,
                          fun           = "enrichGO",
                          keyType = "SYMBOL",
                          pvalueCutoff  = 0.05,
                          pAdjustMethod = "BH",
                          ont = "CC")

zdd_bp <- compareCluster(zdd_genes,
                         OrgDb = org.Mm.eg.db,
                         fun           = "enrichGO",
                         keyType = "SYMBOL",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         ont = "BP")
zdd_mf <- compareCluster(zdd_genes,
                         OrgDb = org.Mm.eg.db,
                         fun           = "enrichGO",
                         keyType = "SYMBOL",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         ont = "MF")
zdd_cc <- compareCluster(zdd_genes,
                         OrgDb = org.Mm.eg.db,
                         fun           = "enrichGO",
                         keyType = "SYMBOL",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         ont = "CC")

zdd_t_bp <- compareCluster(zdd_genes_targets,
                           OrgDb = org.Mm.eg.db,
                           fun           = "enrichGO",
                           keyType = "SYMBOL",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
                           ont = "BP")
zdd_t_mf <- compareCluster(zdd_genes_targets,
                           OrgDb = org.Mm.eg.db,
                           fun           = "enrichGO",
                           keyType = "SYMBOL",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
                           ont = "MF")
zdd_t_cc <- compareCluster(zdd_genes_targets,
                           OrgDb = org.Mm.eg.db,
                           fun           = "enrichGO",
                           keyType = "SYMBOL",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
                           ont = "CC")


# Output results ----------------------------------------------------------
# Npeaks data
write_tsv(invitro_data,"../../results/FinalTables/invitro_table.txt")


# GO data
save(g_bp, g_mf, g_cc,zdd_bp,zdd_mf, zdd_cc,zdd_t_bp,zdd_t_mf, zdd_t_cc,g_cat_bp, g_cat_mf,g_cat_cc,  file = "../../results/invitro/figure_invitro/go_cat.RData")

