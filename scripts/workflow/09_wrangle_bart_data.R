# Author: Melyssa Minto
# This script will visualize BART results of ZicP7 and ZicP60 peaks


# Load libraries ----------------------------------------------------------
library(tidyverse)
library(xlsx)
library(RRHO)

# Read in data ------------------------------------------------------------
result_paths = paste0("../../results/peak_gene/", list.files(path = "../../results/peak_gene/", pattern = "_bart_results.txt", recursive = T))
names(result_paths) = str_extract(result_paths, "(?<=_gene/).+(?=/)")
result_paths2 =  paste0("../../results/DiffExp_ZicChIP/", list.files(path = "../../results/DiffExp_ZicChIP/", pattern = "_bart_results.txt", recursive = T))
names(result_paths2) = paste0(str_extract(result_paths2, "(?<=_ZicChIP/).+(?=/)"), "_all")

all_result_paths = c(result_paths, result_paths2)

bart_data <- all_result_paths %>%
  map_df(read_tsv, .id = "filename")   

geneExp <- read_tsv("../../results/DiffExp_RNA/GeneExp_data.tsv")

tf_gene_mapping = data.frame(TF = c("BIRA", "RPA", "C17orf96", "TP53", rep("TLE", 7), "PR", rep("STAT5",2), "ZNF143", "FAM60A", "ZNF384", "TP63", rep("NP23", 2), "SUPT5H", "NRDC" ),
                             SYMBOL = c("Hlcs", "Rpa1", "Epop","Trp53", paste0("TLE", 1:7), "Pgr", "Stat5a", "Stat5b","Zfp143", "Sinhcaf", "Zfp384","Trp63", "Nup98", "Phf23", "Supt5", "Nrd1"   )
                             
                             )

# Declare functions -------------------------------------------------------


map_bart <-function(cond1, p.val.cut.off){
  bart_data %>% 
    # filter out for user specified conditions
    dplyr::filter(filename %in% cond1) %>% 
    # select the columns needed
    dplyr::select(filename, TF, irwin_hall_pvalue) %>% 
    # pivot data such that each column holds the corresponding p value for condition
    pivot_wider(names_from = filename, values_from = irwin_hall_pvalue) %>%
    # determine if each row is sig in one or both conditions (or not at all)
    dplyr::mutate(bart_sig = ifelse(get(cond1) < p.val.cut.off, "Enriched", "N.S.")) %>% 
    # Rename some of the symbols using tf_gene_mapping
    left_join(tf_gene_mapping) %>% 
    dplyr::mutate(SYMBOL = ifelse(is.na(SYMBOL), str_to_title(TF), str_to_title(SYMBOL))) %>% 
    # split fusion  proteins
    dplyr::mutate(SYMBOL = str_split(as.character(SYMBOL), "(-)(?=[A-Z])")) %>%
    unnest(SYMBOL) %>% 
    # add gene expression data
    left_join(geneExp) 
}

rrho_tf <- function(cond1, cond2, bart_data1, bart_data2, lims = NULL){
  
  # preparing data
  dat1 = bart_data1 %>% 
    dplyr::filter(!grepl("ZIC", TF)) %>% 
    dplyr::select(TF, 2 ) %>% 
    distinct() %>% 
    dplyr::arrange(TF)
  
  colnames(dat1) = c("TF", "FDR")
  
  
  dat2 = bart_data2 %>% 
    dplyr::filter(!grepl("ZIC", TF)) %>% 
    dplyr::select(TF, 2 ) %>% 
    distinct() %>% 
    dplyr::arrange(TF)
  
  
  colnames(dat2) = c("TF", "FDR")
  
  
  
  # making output dir
  system(paste0("rm -r  ../../results/peak_gene/rrho/", cond1,"v" ,cond2))
  system(paste0("mkdir ../../results/peak_gene/rrho/", cond1,"v" ,cond2))
  
  object = RRHO(as.data.frame(dat1), 
                as.data.frame(dat2), 
                labels = c(cond1, cond2), 
                outputdir = paste0("../../results/peak_gene/rrho/",  cond1,"v" ,cond2, "/"), 
                alternative = "enrichment", 
                plots = T,
                BY = T,
                stepsize = 1)
  
  object
}
pull_distinct_enriched <-function(cond1, cond2, bart_data1, bart_data2){
  # get the list of non overlapping tfs 
  overlappingTFs = 
    list.files(paste0("../../results/peak_gene/rrho/", cond1,"v" ,cond2 ), pattern = ".csv", full.names = T) %>% 
    map_df(read_csv, col_names = F) %>% 
    dplyr::filter(!grepl("ZIC", X1)) %>% 
    pull(X1)
  
  allTFs = bart_data1 %>% 
    dplyr::filter(!grepl("ZIC", TF)) %>% 
    dplyr::select(TF) %>% 
    distinct() %>% 
    pull(TF)
  
  non_overlappingTFs = setdiff(allTFs, overlappingTFs)
  
  
  # output the ones that are not enriched 
  distinct_erichedTFs = bind_rows(
    bart_data1 %>% 
      dplyr::filter(TF %in% non_overlappingTFs ) %>% 
      dplyr::filter(get(cond1) < 0.05) %>% 
      dplyr::select("TF", enrich_sig = cond1) %>% 
      distinct() %>% 
      dplyr::mutate(enriched_in = cond1) 
    ,
    bart_data2 %>% 
      dplyr::filter(TF %in% non_overlappingTFs ) %>% 
      dplyr::filter(get(cond2) < 0.05) %>% 
      dplyr::select("TF", enrich_sig = cond2) %>% 
      distinct() %>% 
      dplyr::mutate(enriched_in = cond2)
  ) %>% 
    group_by(TF) %>% 
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

runRRHO <- function(cond1, cond2, bart_data1, bart_data2, lims = NULL){
  # running rrho
  obj = rrho_tf(cond1, cond2, bart_data1, bart_data2, lims)
  # getting distinctly enriched TFs
  distinct_enriched = pull_distinct_enriched(cond1, cond2, bart_data1, bart_data2)
  # rrho heatmap
  pmap = plotRRHO(cond1, cond2, obj, lims)
  
  #making final table
  table = full_join(bart_data1, bart_data2, by = c("TF", "SYMBOL", "p7_mean", "p60_mean", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "gene_sig") ) %>% 
    left_join(distinct_enriched) %>%
    dplyr::mutate(enriched_in = case_when(is.na(enriched_in)  ~ "Similar Enrichment",
                                          !is.na(enriched_in) ~ enriched_in)) %>% 
    dplyr::select("TF","SYMBOL", all_of(cond1), all_of(cond2),"baseMean",  "p60_mean", "p7_mean", "log2FoldChange","gene_sig", "enriched_in") 
  
  list(distinct_enriched = distinct_enriched, pmap = pmap, table = table)
  
  
}


# Wrangle data ------------------------------------------------------------

activating = map_bart("activating", 0.05) %>% write_tsv("../../results/peak_gene/activating/bart_results.txt")
repressive = map_bart("repressive", 0.05) %>% write_tsv("../../results/peak_gene/repressive/bart_results.txt")
late_activating = map_bart("late_activating", 0.05) %>% write_tsv("../../results/peak_gene/late_activating/bart_results.txt")
late_repressive = map_bart("late_repressive", 0.05) %>% write_tsv("../../results/peak_gene/late_repressive/bart_results.txt")
early_activating = map_bart("early_activating", 0.05) %>% write_tsv("../../results/peak_gene/early_activating/bart_results.txt")
early_repressive = map_bart("early_repressive", 0.05) %>% write_tsv("../../results/peak_gene/early_repressive/bart_results.txt")
early = map_bart("early", 0.05) %>% write_tsv("../../results/peak_gene/early/bart_results.txt")
late = map_bart("late", 0.05) %>% write_tsv("../../results/peak_gene/late/bart_results.txt")
early_all = map_bart("early_all", 0.05) %>% write_tsv("../../results/DiffExp_ZicChIP/early/bart_results.txt")
late_all = map_bart("late_all", 0.05) %>% write_tsv("../../results/DiffExp_ZicChIP/late/bart_results.txt")
static_all = map_bart("late_all", 0.05) %>% write_tsv("../../results/DiffExp_ZicChIP/static/bart_results.txt")

earlyVlate_act = runRRHO("early_activating", "late_activating", early_activating, late_activating, lims = c(0, 400))
early_actVrep = runRRHO("early_activating", "early_repressive", early_activating, early_repressive, lims = c(0, 400))
earlyVlate_rep = runRRHO("early_repressive", "late_repressive", early_repressive,late_repressive, lims = c(0, 400))
late_actvrep = runRRHO("late_activating", "late_repressive", late_activating, late_repressive, lims = c(0, 400))
actvrep = runRRHO("activating", "repressive", activating, repressive, lims = c(0, 400))
earlyvlate = runRRHO("early", "late", early, late, lims = c(0, 400))
early_allVloop = runRRHO("early_all", "early", early_all, early, lims = c(0, 400))
late_allVloop = runRRHO("late_all", "late", late_all, late, lims = c(0, 400))
earlyVlate_all = runRRHO("early_all", "late_all", early_all, late_all, lims = c(0,400))
# Save table --------------------------------------------------------------

earlyVlate_act$table %>% write_tsv("../../results/peak_gene/rrho/early_late_activating.txt")
early_actVrep$table %>% write_tsv("../../results/peak_gene/rrho/early_activating_repressive.txt")
earlyVlate_rep$table %>% write_tsv("../../results/peak_gene/rrho/early_late_repressive.txt")
late_actvrep$table %>% write_tsv("../../results/peak_gene/rrho/late_activating_repressive.txt")
actvrep$table %>% write_tsv("../../results/peak_gene/rrho/activating_repressive.txt")
earlyvlate$table %>% write_tsv("../../results/peak_gene/rrho/early_late.txt")
early_allVloop$table %>% write_tsv("../../results/peak_gene/rrho/early_allvloop.txt")
late_allVloop$table %>% write_tsv("../../results/peak_gene/rrho/late_allvloop.txt")
earlyVlate_all$table %>% write_tsv("../../results/peak_gene/rrho/earlyvlate_all.txt")
