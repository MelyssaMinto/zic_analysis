# Author: Melyssa Minto
# This script will compare homer results of ZicP7 and ZicP60 peaks

# Load libraries ----------------------------------------------------------
library(tidyverse)
library(xlsx)
library(RRHO)

# Read in data ------------------------------------------------------------

homer_files = list.files("../../results/homer_results/", pattern = "homer_results.tsv", recursive = T) 
name = str_extract(homer_files, "[^/]+")
homer_files = paste0("../../results/homer_results/", homer_files) 
names(homer_files) = name 


homer_data <- homer_files %>%
  map_df(read_tsv, .id = "filename") %>% 
  dplyr::mutate(filename = case_when( filename == "P60_peaks_DOWNGenes" ~ "late_repressive",
                                      filename == "P60_peaks_UpGenes" ~ "late_activating",
                                      filename == "P7_peaks_DOWNGenes" ~ "early_activating",
                                      filename == "P7_peaks_UpGenes" ~ "early_repressive",
                                      filename == "P60vP7_DOWN" ~ "early_all",
                                      filename == "P60vP7_UP" ~ "late_all",
                                      filename == "P60vP7_NS" ~ "static_all",
                                      filename %in% c("early", "activating", "repressive", "late") ~ filename
                                      
                                      ))

# Functions ---------------------------------------------------------------

rrho_tf <- function(cond1, cond2, homer_data1, homer_data2, lims = NULL){
  
  # preparing data
  dat1 = homer_data1 %>% 
    # removing enriched motifs in the zic family
    dplyr::filter(!grepl("ZIC", `Motif Name`, ignore.case = T)) %>% 
    # getting the enriched TFs and their FDR
    dplyr::select(`Motif Name`, `q-value (Benjamini)`  ) %>% 
    distinct() %>% 
    # if there are duplicated in enriched TF,take the max FDR
    group_by(`Motif Name`) %>% 
    dplyr::slice_max( `q-value (Benjamini)`) %>% 
    ungroup()
  colnames(dat1) = c("TF", "FDR")
  
  
  dat2 = homer_data2 %>% 
    # removing enriched motifs in the zic family
    dplyr::filter(!grepl("ZIC", `Motif Name`, ignore.case = T)) %>% 
    # getting the enriched TFs and their FDR
    dplyr::select(`Motif Name`, `q-value (Benjamini)`  ) %>% 
    distinct() %>% 
    # if there are duplicated in enriched TF,take the max FDR
    group_by(`Motif Name`) %>% 
    dplyr::slice_max( `q-value (Benjamini)`) %>% 
    ungroup
  colnames(dat2) = c("TF", "FDR")
  
  
  # making output dir
  system(paste0("rm -r  ../../results/peak_gene/rrho/", cond1,"v" ,cond2, "_homer"))
  system(paste0("mkdir ../../results/peak_gene/rrho/", cond1,"v" ,cond2, "_homer"))
  
  object = RRHO(as.data.frame(dat1), 
                as.data.frame(dat2), 
                labels = c(cond1, cond2), 
                outputdir = paste0("../../results/peak_gene/rrho/",  cond1,"v" ,cond2, "_homer/"), 
                alternative = "enrichment", 
                plots = T,
                BY = T,
                stepsize = 1)
  
  object
}
pull_distinct_enriched <-function(cond1, cond2, homer_data1, homer_data2){
  # get the list of non overlapping tfs 
  overlappingTFs = 
    list.files(paste0("../../results/peak_gene/rrho/", cond1,"v" ,cond2, "_homer" ), pattern = ".csv", full.names = T) %>% 
    map_df(read_csv, col_names = F) %>% 
    dplyr::filter(!grepl("ZIC", X1)) %>% 
    pull(X1)
  
  allTFs = homer_data1 %>% 
    dplyr::filter(!grepl("ZIC", `Motif Name`, ignore.case = T)) %>% 
    dplyr::select(`Motif Name`) %>% 
    distinct() %>% 
    pull(`Motif Name`)
  
  non_overlappingTFs = setdiff(allTFs, overlappingTFs)
  
  
  # output the ones that are not enriched 
  distinct_erichedTFs = bind_rows(
    homer_data1 %>% 
      dplyr::filter(`Motif Name` %in% non_overlappingTFs ) %>% 
      dplyr::filter(`q-value (Benjamini)`< 0.05) %>% 
      dplyr::select("Motif Name", enrich_sig = `q-value (Benjamini)`) %>% 
      distinct() %>% 
      dplyr::mutate(enriched_in = cond1) 
    ,
    homer_data2 %>% 
      dplyr::filter(`Motif Name` %in% non_overlappingTFs ) %>% 
      dplyr::filter(`q-value (Benjamini)` < 0.05) %>% 
      dplyr::select("Motif Name", enrich_sig = `q-value (Benjamini)`) %>% 
      distinct() %>% 
      dplyr::mutate(enriched_in = cond2)
  ) %>% 
    group_by(`Motif Name`) %>% 
    dplyr::slice_min(enrich_sig) %>% 
    ungroup()
  
  
}

plotRRHO <-function(cond1, cond2, object, lims){
  
  # extracting and rotating rrho matrix
  mat = object$hypermat
  mat = t(mat)
  
  # setting p values that are Inf to hightest p value (countable)
  if(max(mat) == Inf){
    cat("\nWARNING: some p-values are Inf, so they are being replaced with the max on the color bar\n")
    mat[which(mat %in% Inf)] = max( mat[mat!=max(mat)] )
    
  }
  # plotting RRHO
  mycol <- c("navy", "blue", "cyan", "green", "yellow", "red", "red4")
  
  pmap = as.data.frame(mat) %>% 
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
    labs(x = cond1, y = cond2, fill = "-log10(p-value)") +
    ggtitle("Rank Rank Hypergeometric Overlap Map")+
    theme(plot.title = element_text(hjust = 0.5, face="bold"))
  
  pmap
}

runRRHO <- function(cond1, cond2, lims = NULL){
  
  data1 = homer_data %>% dplyr::filter(filename == cond1)
  data2 = homer_data %>% dplyr::filter(filename == cond2)
   
  # running rrho
  obj = rrho_tf(cond1, cond2, data1 , data2, lims)
  # getting distinctly enriched TFs
  distinct_enriched = pull_distinct_enriched(cond1, cond2, data1, data2)
  # rrho heatmap
  pmap = plotRRHO(cond1, cond2, obj, lims)
  
  #making final table
  table = full_join(data1, data2, by = c("Motif Name", "SYMBOL","q-value (Benjamini)", "p7_mean", "p60_mean", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "gene_sig") ) %>% 
    left_join(distinct_enriched) %>%
    dplyr::mutate(enriched_in = case_when(is.na(enriched_in)  ~ "Similar Enrichment",
                                          !is.na(enriched_in) ~ enriched_in)) %>% 
    dplyr::select("Motif Name","SYMBOL", "q-value (Benjamini)","baseMean",  "p60_mean", "p7_mean", "log2FoldChange","gene_sig", "enriched_in") 
  
  list(distinct_enriched = distinct_enriched, pmap = pmap, table = table)
  
  
}


# Wrangle data ------------------------------------------------------------

earlyVlate_act = runRRHO(cond1 = "early_activating", cond2 = "late_activating", lims = NULL)
early_actVrep = runRRHO(cond1 = "early_activating", cond2 = "early_repressive", lims = NULL)
earlyVlate_rep = runRRHO(cond1 = "early_repressive", cond2 = "late_repressive", lims = NULL)
late_actvrep = runRRHO(cond1 = "late_activating", cond2 = "late_repressive", lims = NULL)
actvrep = runRRHO(cond1 = "activating", cond2 = "repressive", lims = NULL)
earlyvlate = runRRHO(cond1 = "early", cond2 = "late", lims = NULL)
early_allVloop = runRRHO(cond1 = "early_all", cond2 = "early", lims = NULL)
late_allVloop = runRRHO(cond1 = "late_all", cond2 = "late", lims = NULL)
earlyVlate_all = runRRHO(cond1 = "early_all", cond2 = "late_all", lims = NULL)


# Save table --------------------------------------------------------------

earlyVlate_act$table %>% write_tsv("../../results/peak_gene/rrho/early_late_activating_homer.txt")
early_actVrep$table %>% write_tsv("../../results/peak_gene/rrho/early_activating_repressive_homer.txt")
earlyVlate_rep$table %>% write_tsv("../../results/peak_gene/rrho/early_late_repressive_homer.txt")
late_actvrep$table %>% write_tsv("../../results/peak_gene/rrho/late_activating_repressive_homer.txt")
actvrep$table %>% write_tsv("../../results/peak_gene/rrho/activating_repressive_homer.txt")
earlyvlate$table %>% write_tsv("../../results/peak_gene/rrho/early_late_homer.txt")
early_allVloop$table %>% write_tsv("../../results/peak_gene/rrho/early_allvloop_homer.txt")
late_allVloop$table %>% write_tsv("../../results/peak_gene/rrho/late_allvloop_homer.txt")
earlyVlate_all$table %>% write_tsv("../../results/peak_gene/rrho/earlyvlate_all_homer.txt")
