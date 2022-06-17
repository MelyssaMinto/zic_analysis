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



# functions ---------------------------------------------------------------

customRRHO <- function(data, dir, conditions, label, lims=NULL){
  
  data = data %>% 
    dplyr::select(matches(c("names","zic[1|2]_log2FoldChange","zic[1|2]_padj","dev_log2FoldChange","dev_padj"  )) ) %>% 
    distinct() %>% 
    drop_na()
    
  zic_signal = data %>% 
    dplyr::select(matches(c("names","zic[1|2]_log2FoldChange","zic[1|2]_padj"  )) ) %>% 
    dplyr::mutate(signal = sign(.[[2]]) * log10(.[[3]]) ) %>% 
    dplyr::select(ends_with("names"), signal) %>% 
  distinct() 
  
  dev_signal = data %>% 
    dplyr::select(matches(c("names", "dev_log2FoldChange","dev_padj"  )) ) %>% 
    dplyr::mutate(signal = sign(.[[2]]) * log10(.[[3]]) ) %>% 
    dplyr::select(ends_with("names"), signal) %>% 
    distinct() 
  
  dev_signal_neg = dev_signal %>% 
    dplyr::mutate(signal = -signal )
  
  object = RRHO(as.data.frame(zic_signal),
                       as.data.frame(dev_signal),
                       BY=TRUE,
                       stepsize = 100,
                       alternative='two.sided',
                       log10.ind = TRUE,
                       plot = TRUE,
                       outputdir = "../../results/invitro/rrho/",
                       labels =  c(label,"D7_D3"))
  
  object_neg = RRHO(as.data.frame(zic_signal),
                as.data.frame(dev_signal_neg),
                BY=TRUE,
                stepsize = 100,
                alternative='two.sided',
                log10.ind = TRUE,
                plot = TRUE,
                outputdir = "../../results/invitro/rrho/",
                labels = c(label,"D3_D7"))
  
  
  # extracting and rotating matrix
  mat = object$hypermat
  mat_neg = object_neg$hypermat
  
  mat = t(mat)
  mat_neg = t(mat_neg)
  
  # setting p values that are Inf to hightest p value (countable)
  if(max(mat) == Inf){
    cat("\nWARNING: some p-values are Inf, so they are being replaced with the max on the color bar\n")
    mat[which(mat %in% Inf)] = max( mat[mat!=max(mat)] )
    
  }
  
  if(max(mat_neg) == Inf){
    cat("\nWARNING: some p-values are Inf, so they are being replaced with the max on the color bar\n")
    mat_neg[which(mat_neg %in% Inf)] = max( mat_neg[mat_neg!=max(mat_neg)] )
    
  }
  
  # adding the signs
  mat_signs = t(object$hypermat.signs)
  mat = as.data.frame(mat)*as.data.frame(mat_signs)
  
  mat_neg_signs = t(object_neg$hypermat.signs)
  mat_neg = as.data.frame(mat_neg)*as.data.frame( mat_neg_signs)
  
  # plot using ggplot
  mycol <- c("navy", "blue", "cyan", "green", "yellow", "red", "red4")
  
  p_map = mat %>%
    rownames_to_column("diff_sig1") %>%
    pivot_longer(-c(diff_sig1), names_to = "diff_sig2", values_to = "hypergeometric_pvals") %>% 
    dplyr::mutate(diff_sig1= fct_relevel(diff_sig1,rownames(as.data.frame(mat)))) %>%
    dplyr::mutate(diff_sig2= fct_relevel(diff_sig2,colnames(as.data.frame(mat))))%>%
    ggplot(aes(x=diff_sig2, y=diff_sig1, fill=hypergeometric_pvals)) + 
    geom_raster() +
    scale_fill_gradientn(colours = mycol, limits = lims, n.breaks=7)+
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          legend.key.height= unit(2.7, 'cm'),
          legend.key.width= unit(1, 'cm')) +
    labs(x = "", y = "", fill = "-log10(p-value)") +
    ggtitle("Rank Rank Hypergeometric Overlap Map")+
    theme(plot.title = element_text(hjust = 0.5, face="bold"))
  
  p_map_neg = mat_neg %>% 
    rownames_to_column("diff_sig1") %>%
    pivot_longer(-c(diff_sig1), names_to = "diff_sig2", values_to = "hypergeometric_pvals") %>% 
    dplyr::mutate(diff_sig1= fct_relevel(diff_sig1,rownames(as.data.frame(mat_neg)))) %>%
    dplyr::mutate(diff_sig2= fct_relevel(diff_sig2,colnames(as.data.frame(mat_neg))))%>%
    ggplot(aes(x=diff_sig2, y=diff_sig1, fill=hypergeometric_pvals)) + 
    geom_raster() +
    scale_fill_gradientn(colours = mycol, limits = lims, n.breaks=7)+
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          legend.key.height= unit(2.7, 'cm'),
          legend.key.width= unit(1, 'cm')) +
    labs(x = "", y = "", fill = "-log10(p-value)") +
    ggtitle("Rank Rank Hypergeometric Overlap Map")+
    theme(plot.title = element_text(hjust = 0.5, face="bold"))
  
  
  
  
  png(paste0(dir, "RRHO.png"), res = 300, units = 'in', height = 6.62,width =7.02 )
  print( p_map  )
  
  grid.lines(x = unit(c(0.05,.8), "npc"), y = unit(c(0.03,0.03), "npc"), arrow = arrow(length = unit(.4, "cm"), type = "closed", ends = "both"))
  grid.text(conditions[2], x = unit(0.7, "npc"), y = unit(0.015, "npc"))
  grid.text(conditions[1], x = unit(0.15, "npc"), y = unit(0.015, "npc"))
  
  grid.lines(x = unit(c(0.03,.03), "npc"), y = unit(c(0.06,.94), "npc"), arrow = arrow(length = unit(.4, "cm"), type = "closed", ends = "both"), name = "Help")
  grid.text(conditions[4], x = unit(0.015, "npc"), y = unit(0.85, "npc"),rot = 90)
  grid.text(conditions[3], x = unit(0.015, "npc"), y = unit(0.17, "npc"),rot = 90)
  
  dev.off()
  
  png(paste0(dir, "RRHO_neg.png"), res = 300, units = 'in', height = 6.62,width =7.02 )
  print( p_map_neg  )
  
  grid.lines(x = unit(c(0.05,.8), "npc"), y = unit(c(0.03,0.03), "npc"), arrow = arrow(length = unit(.4, "cm"), type = "closed", ends = "both"))
  grid.text(conditions[2], x = unit(0.7, "npc"), y = unit(0.015, "npc"))
  grid.text(conditions[1], x = unit(0.15, "npc"), y = unit(0.015, "npc"))
  
  grid.lines(x = unit(c(0.03,.03), "npc"), y = unit(c(0.06,.94), "npc"), arrow = arrow(length = unit(.4, "cm"), type = "closed", ends = "both"), name = "Help")
  grid.text(conditions[3], x = unit(0.015, "npc"), y = unit(0.85, "npc"),rot = 90)
  grid.text(conditions[4], x = unit(0.015, "npc"), y = unit(0.17, "npc"),rot = 90)
  
  dev.off()
}

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




# Run RRHO analysis -------------------------------------------------------

customRRHO(knockdown_genes_zic1_dev,  "../../results/invitro/rrho/zic1_dev_", conditions = c("Zic1 KD - Up", "Zic1 KD - Down", "D7", "D3"), label = "Zic1KD v WT", lims = NULL )
customRRHO(knockdown_genes_zic2_dev,  "../../results/invitro/rrho/zic2_dev_", conditions = c("Zic2 KD - Up", "Zic1 KD - Down", "D7", "D3"), label = "Zic2KD v WT", lims = NULL)
