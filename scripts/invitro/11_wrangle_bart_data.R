# Author: Melyssa Minto
# This script will visualize BART results of ZicP7 and ZicP60 peaks


# Load libraries ----------------------------------------------------------
library(tidyverse)
library(xlsx)
library(RRHO)

# Read in data ------------------------------------------------------------
result_paths = paste0("../../results/invitro/bart_results/", list.files(path = "../../results/invitro/bart_results/", pattern = "_bart_results.txt", recursive = T))
names(result_paths) = str_extract(result_paths, "(?<=_results/).+(?=/)")

bart_data <- result_paths %>%
  map_df(read_tsv, .id = "filename")   

geneExp <- read_tsv("../../results/invitro/diffexpr/res_d7_d3_ashr.tsv") %>% 
  dplyr::rename(SYMBOL = gene_names)

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
  system(paste0("rm -r  ../../results/invitro/beds/rrho/", cond1,"v" ,cond2))
  system(paste0("mkdir ../../results/invitro/beds/rrho/", cond1,"v" ,cond2))
  
  object = RRHO(as.data.frame(dat1), 
                as.data.frame(dat2), 
                labels = c(cond1, cond2), 
                outputdir = paste0("../../results/invitro/beds/rrho/",  cond1,"v" ,cond2, "/"), 
                alternative = "enrichment", 
                plots = T,
                BY = T,
                stepsize = 1)
  
  object
}
pull_distinct_enriched <-function(cond1, cond2, bart_data1, bart_data2){
  # get the list of non overlapping tfs 
  overlappingTFs = 
    list.files(paste0("../../results/invitro/beds/rrho/", cond1,"v" ,cond2 ), pattern = ".csv", full.names = T) %>% 
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

plotRRHO <-function(cond1, cond2, name1 = NULL, name2=NULL, object, lims){
  if(is.null(name1) | is.null(name2)){
    name1 = cond1
    name2 = cond2
  }
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
          legend.key.height= unit(1, 'cm'),
          legend.key.width= unit(.5, 'cm')) +
    labs(x = name1, y = name2, fill = "-log10(p-value)") +
    theme(plot.title = element_text(hjust = 0.5, face="bold"))
  
  pmap
}

runRRHO <- function(cond1, cond2, name1 = NULL, name2 = NULL, bart_data1, bart_data2, lims = NULL){
  # running rrho
  obj = rrho_tf(cond1, cond2, bart_data1, bart_data2, lims)
  # getting distinctly enriched TFs
  distinct_enriched = pull_distinct_enriched(cond1, cond2, bart_data1, bart_data2)
  # rrho heatmap
  pmap = plotRRHO(cond1, cond2, name1, name2, obj, lims)
  
  #making final table
  table = full_join(bart_data1, bart_data2, by = c("TF", "SYMBOL",  "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj") ) %>% 
    left_join(distinct_enriched) %>%
    dplyr::mutate(enriched_in = case_when(is.na(enriched_in)  ~ "Similar Enrichment",
                                          !is.na(enriched_in) ~ enriched_in)) %>% 
    dplyr::select("TF","SYMBOL", all_of(cond1), all_of(cond2),"baseMean", "log2FoldChange", "enriched_in") 
  
  list(distinct_enriched = distinct_enriched, pmap = pmap, table = table)
  
  
}


# Wrangle data ------------------------------------------------------------
dev_anchors = map_bart("dev_anchors", 0.05) %>% write_tsv("../../results/invitro/bart_results/dev_anchors/bart_results.txt")
dev_zicPeaks = map_bart("dev_zicPeaks", 0.05) %>% write_tsv("../../results/invitro/bart_results/dev_zicPeaks/bart_results.txt")
zd_anchors = map_bart("zd_anchors", 0.05) %>% write_tsv("../../results/invitro/bart_results/zd_anchors/bart_results.txt")
zd_zicPeaks = map_bart("zd_zicPeaks", 0.05) %>% write_tsv("../../results/invitro/bart_results/zd_zicPeaks/bart_results.txt")
zdd_anchors = map_bart("zdd_anchors", 0.05) %>% write_tsv("../../results/invitro/bart_results/zdd_anchors/bart_results.txt")
zdd_zicPeaks = map_bart("zdd_zicPeaks", 0.05) %>% write_tsv("../../results/invitro/bart_results/zdd_zicPeaks/bart_results.txt")


devVzd_anchors = runRRHO("dev_anchors", "zd_anchors", "Developmental", "Zic Dependent", dev_anchors, zd_anchors, lims = c(0, 340))
devVzdd_anchors = runRRHO("dev_anchors", "zdd_anchors",  "Developmental", "Zic Dependent Developmental",dev_anchors, zdd_anchors, lims = c(0, 340))
zdVzdd_anchors = runRRHO("zd_anchors", "zdd_anchors",  "Zic Dependent", "Zic Dependent Developmental",zd_anchors, zdd_anchors, lims = c(0, 340))

devVzd_zicPeaks = runRRHO("dev_zicPeaks", "zd_zicPeaks", "Developmental", "Zic Dependent", dev_zicPeaks, zd_zicPeaks, lims = c(0, 340))
devVzdd_zicPeaks = runRRHO("dev_zicPeaks", "zdd_zicPeaks",  "Developmental", "Zic Dependent",dev_zicPeaks, zdd_zicPeaks, lims = c(0, 340))
zdVzdd_zicPeaks = runRRHO("zd_zicPeaks", "zdd_zicPeaks",  "Zic Dependent", "Zic Dependent Developmental",zd_zicPeaks, zdd_zicPeaks, lims = c(0, 340))



# Save table --------------------------------------------------------------

devVzd_anchors$table %>% write_tsv("../../results/invitro/beds/rrho/devVzd_anchors.txt")
devVzdd_anchors$table %>% write_tsv("../../results/invitro/beds/rrho/devVzdd_anchors.txt")
zdVzdd_anchors$table %>% write_tsv("../../results/invitro/beds/rrho/zdVzdd_anchors.txt")

devVzd_zicPeaks$table %>% write_tsv("../../results/invitro/beds/rrho/devVzd_zicPeaks.txt")
devVzdd_zicPeaks$table %>% write_tsv("../../results/invitro/beds/rrho/devVzdd_zicPeaks.txt")
zdVzdd_zicPeaks$table %>% write_tsv("../../results/invitro/beds/rrho/zdVzdd_zicPeaks.txt")

# Save heatmaps -----------------------------------------------------------

png("../../figures/devVzd_anchors.png", units = "in", height = 3, width = 4, res = 300)
devVzd_anchors$pmap 
dev.off()

png("../../figures/devVzdd_anchors.png", units = "in", height = 3, width = 4, res = 300)
devVzdd_anchors$pmap
dev.off()


png("../../figures/zdVzdd_anchors.png", units = "in", height = 3, width = 4, res = 300)
zdVzdd_anchors$pmap
dev.off()

png("../../figures/devVzd_zicPeaks.png", units = "in", height = 3, width = 4, res = 300)
devVzd_zicPeaks$pmap 
dev.off()

png("../../figures/devVzdd_zicPeaks.png", units = "in", height = 3, width = 4, res = 300)
devVzdd_zicPeaks$pmap
dev.off()


png("../../figures/zdVzdd_zicPeaks.png", units = "in", height = 3, width = 4, res = 300)
zdVzdd_zicPeaks$pmap
dev.off()


