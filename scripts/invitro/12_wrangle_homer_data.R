# Author: Melyssa Minto
# This script will compare homer results of ZicP7 and ZicP60 peaks

# Load libraries ----------------------------------------------------------
library(RRHO)
library(tidyverse)
library(readxl)
library(reshape2)
library(foreach)
library(xlsx)

# Read in data ------------------------------------------------------------

# peak files
peakFilesOfInterest = gsub(".bed","" , list.files("../../results/invitro/beds/", pattern = ".bed", recursive = T))[c(1,3,5)]

# motif enrichment results
homerResultsPath = paste0("../../results/invitro/homer_results/",peakFilesOfInterest, "/results/invitro/knownResults.txt")
homerResults <- homerResultsPath %>%
  map(read_tsv) 


# DiffExp data
DEResults <- read_tsv("../../results/invitro/diffexpr/res_d7_d3_ashr.tsv") %>% 
  dplyr::rename(SYMBOL = gene_names) %>% 
  dplyr::mutate(gene_sig = case_when(log2FoldChange > 1 & padj < 0.05 ~ "UP",
                log2FoldChange < -1 & padj < 0.05 ~ "DOWN",
                is.na(padj) ~ "filtered",
                TRUE ~ "N.S."))

# protein2gene name conversion
proteinGene = read_csv("../../../../../genomeData/mm10/proteinsThatAreNotNamedAsTheirGeneNamesWithGeneNames.csv", col_names = F)

# Define Functions --------------------------------------------------------

clean_name <- function(df, gene_exp){
  # this function parses through the protein names and extract the name and family
  # input
  #   df =data.frame of homer outputs
  #   gene.exp = deseq results with gene name
  # output
  #   cleaned homer output with Tf name and TF family 
  
  # wrangle
  gene_exp = gene_exp %>% 
    dplyr::mutate(SYMBOL = tolower(SYMBOL))
  
  # extracting name and family from the homer output
  name_fam= sub("\\/.*", "", df$`Motif Name`)
  name=sub("\\(.*", "",name_fam)
  fam=sub(".*\\((.*)\\).*", "\\1", name_fam, perl=TRUE) 
  fam[!grepl("\\(", name_fam)] <- NA
  
  df$name=name
  df$fam=fam
  
  #creating the gene names column for the corrsponding protien name
  df$SYMBOL = tolower(df$name)
  df$SYMBOL = gsub("znf", "zfp", df$SYMBOL, ignore.case = T)  # all zinc fingers are prfixed with zfp
  df$SYMBOL = gsub("\\.", "-", df$SYMBOL, ignore.case = T)  
  df$SYMBOL = gsub("-fusion", "", df$SYMBOL, ignore.case = T)  
  df$SYMBOL = gsub("-distal", "", df$SYMBOL, ignore.case = T)  
  df$SYMBOL = gsub("-halfsite", "", df$SYMBOL, ignore.case = T)  
  df$SYMBOL = gsub(":ebox", "", df$SYMBOL, ignore.case = T)  
  df$SYMBOL = gsub(":e-box", "", df$SYMBOL, ignore.case = T)  
  df$SYMBOL = gsub("-satelliteelement", "", df$SYMBOL, ignore.case = T)  
  df$SYMBOL = gsub("-short", "", df$SYMBOL, ignore.case = T)  
  df$SYMBOL = gsub("\\+1bp", "", df$SYMBOL, ignore.case = T)  
  
  
  # splitting protein names if they are fused protiens (ie. AP1-Jun will result in 2 rows, one for AP1 and one for jun)
  df_clean = foreach(i = 1:nrow(df), .combine = "rbind") %do%
    {
      
      
      elem = df[i,]
      
      # split ":", "+", " "
      if(grepl(":", elem$SYMBOL)){
        split_elem = rbind(elem, elem)
        split_elem$SYMBOL[1] = sub(":.*", "", elem$SYMBOL)
        split_elem$SYMBOL[2] = sub(".*:", "", elem$SYMBOL)
        new_elem = split_elem
      }else if(grepl("\\+", elem$SYMBOL)){
        split_elem = rbind(elem, elem)
        split_elem$SYMBOL[1] = sub("\\+.*", "", elem$SYMBOL)
        split_elem$SYMBOL[2] = sub(".*\\+", "", elem$SYMBOL)
        new_elem = split_elem
      }else if(grepl(" ", elem$SYMBOL)){
        split_elem = rbind(elem, elem)
        split_elem$SYMBOL[1] = sub(" .*", "", elem$SYMBOL)
        split_elem$SYMBOL[2] = sub(".* ", "", elem$SYMBOL)
        new_elem = split_elem
      }else{
        new_elem = elem
      }
      
      # matching genes to protein
      #   final_elem = foreach(j= 1:nrow(new_elem), .combine= 'rbind') %do%
      #     {
      #       return(left_join(new_elem, gene_exp))
      #       
      #     }
      #   
      
      # subbing gene names for protein names that are incongruent 
      final_elem = foreach( j = 1:nrow(new_elem), .combine = "rbind") %do%
        {
          
          
          
          # extract the gene(s) that correspond to that protein
          genes = ifelse(is_null(protein2gene[[ new_elem$SYMBOL[j] ]]),  new_elem$SYMBOL, protein2gene[[ new_elem$SYMBOL[j] ]])
          
          # repeat the rows for as many genes there are
          fin_elem = new_elem[j,] %>% dplyr::slice(rep(1:n(), each = length(genes)))
          
          # replace the protein name with the extracted gene name(s)
          fin_elem$SYMBOL = genes
          #cat('i=',i, ' j=' , j, '\n')
          return(left_join(fin_elem, gene_exp, by=  "SYMBOL", copy=T))
        }
      
      
      
      return(final_elem)
    }
  
  return(df_clean)
}


# Data Cleaning -----------------------------------------------------------
# > RNA data --------------------------------------------------------------


# annotating filtered genes in significance columns
expData =DEResults %>% 
  dplyr::mutate(sig = ifelse(is.na(padj), "filtered", gene_sig))

# > making proteinGene a dictionary ----------------------------------------

protein2gene = list() # make empty list

# first make a mapping of protein to gene names assuming that the protein and gene symbols are the smae
for( i in 1:nrow(expData)){ 
  protein2gene[[ tolower(expData$SYMBOL[i]) ]] = tolower(expData$SYMBOL[i])
}

# go through the list of protein names that are incomgruent and add the correct gene name(s)
for( i in 1:nrow(proteinGene)){
  prot = trimws(tolower(as.character(proteinGene[i,1])))
  gene = tolower(as.character(proteinGene[i,2]))
  genes = trimws(unlist(strsplit(gene, "\\|")))
  
  if(length(genes) > 1){
    
    protein2gene[[prot]] = genes
  }else{
    protein2gene[[prot]] = gene
  }
  
}

length(unique(names(protein2gene)))

# > Homer data -------------------------------------------------------

# parsing out TF name
homerResultsClean = list()
for(i in 1:length(homerResults)){
  homerResultsClean[[i]]= clean_name(homerResults[[i]], expData)
}

names(homerResultsClean) =  gsub(".bed", "", peakFilesOfInterest) 





# Write data --------------------------------------------------------------


write_tsv(as.data.frame(homerResultsClean[[1]]),paste0("../../results/invitro/homer_results/", names(homerResultsClean)[1], "/homer_results.tsv"))
write_tsv(as.data.frame(homerResultsClean[[2]]),paste0("../../results/invitro/homer_results/", names(homerResultsClean)[2], "/homer_results.tsv"))
write_tsv(as.data.frame(homerResultsClean[[3]]),paste0("../../results/invitro/homer_results/", names(homerResultsClean)[3], "/homer_results.tsv"))


rm(ls())

# Read in data ------------------------------------------------------------

homer_files = list.files("../../results/invitro/homer_results/", pattern = "homer_results.tsv", recursive = T) 
name = str_extract(homer_files, "[^/]+")
homer_files = paste0("../../results/invitro/homer_results/", homer_files) 
names(homer_files) = name 


homer_data <- homer_files %>%
  map_df(read_tsv, .id = "filename") 

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
  system(paste0("rm -r  ../../results/invitro/beds/rrho/", cond1,"v" ,cond2, "_homer"))
  system(paste0("mkdir ../../results/invitro/beds/rrho/", cond1,"v" ,cond2, "_homer"))
  
  object = RRHO(as.data.frame(dat1), 
                as.data.frame(dat2), 
                labels = c(cond1, cond2), 
                outputdir = paste0("../../results/invitro/beds/rrho/",  cond1,"v" ,cond2, "_homer/"), 
                alternative = "enrichment", 
                plots = T,
                BY = T,
                stepsize = 1)
  
  object
}
pull_distinct_enriched <-function(cond1, cond2, homer_data1, homer_data2){
  # get the list of non overlapping tfs 
  overlappingTFs = 
    list.files(paste0("../../results/invitro/beds/rrho/", cond1,"v" ,cond2, "_homer" ), pattern = ".csv", full.names = T) %>% 
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

runRRHO <- function(cond1, cond2, name1=NULL, name2=NULL, lims = NULL){
  
  data1 = homer_data %>% dplyr::filter(filename == cond1)
  data2 = homer_data %>% dplyr::filter(filename == cond2)
   
  # running rrho
  obj = rrho_tf(cond1, cond2, data1 , data2, lims)
  # getting distinctly enriched TFs
  distinct_enriched = pull_distinct_enriched(cond1, cond2, data1, data2)
  # rrho heatmap
  pmap = plotRRHO(cond1, cond2, name1, name2, obj, lims)
  
  #making final table
  table = homer_data %>%
    dplyr::filter(filename %in% c(cond1, cond2)) %>% 
    dplyr::select("Motif Name", SYMBOL, name, filename, "q-value (Benjamini)") %>% 
    distinct() %>% 
    group_by(`Motif Name`, filename) %>% 
    slice_max(`q-value (Benjamini)`) %>% 
    ungroup() %>% 
    
    pivot_wider(id_cols = c(`Motif Name`, SYMBOL, name), names_from = filename, values_from = "q-value (Benjamini)") %>% 
    distinct() %>% 
    left_join(homer_data %>% dplyr::select(SYMBOL, baseMean, log2FoldChange) %>% distinct()) %>% 
    left_join(distinct_enriched) %>% 
    dplyr::mutate(enriched_in = case_when(is.na(enriched_in)  ~ "Similar Enrichment",
                                          !is.na(enriched_in) ~ enriched_in)) %>% 
    dplyr::mutate(SYMBOL = str_to_title(SYMBOL)) %>% 
    dplyr::rename(TF = name)

  return(list(distinct_enriched = distinct_enriched, pmap = pmap, table = table))
  
  
}



# Wrangle data ------------------------------------------------------------

devVzd_anchors = runRRHO("dev_anchors", "zd_anchors", "Developmental", "Zic Dependent", lims = c(0, 15))
devVzdd_anchors = runRRHO("dev_anchors", "zdd_anchors",  "Developmental", "Zic Dependent Developmental",lims = c(0, 15))
zdVzdd_anchors = runRRHO("zd_anchors", "zdd_anchors",  "Zic Dependent", "Zic Dependent Developmental",lims = c(0, 340))


# Save table --------------------------------------------------------------

devVzd_anchors$table %>% write_tsv("../../results/invitro/beds/rrho/devVzd_anchors_homer.txt")
devVzdd_anchors$table %>% write_tsv("../../results/invitro/beds/rrho/devVzdd_anchors_homer.txt")
zdVzdd_anchors$table %>% write_tsv("../../results/invitro/beds/rrho/zdVzdd_anchors_homer.txt")



# Save heatmaps -----------------------------------------------------------

png("../../figures/devVzd_anchors_homer.png", units = "in", height = 3, width = 4, res = 300)
devVzd_anchors$pmap 
dev.off()

png("../../figures/devVzdd_anchors_homer.png", units = "in", height = 3, width = 4, res = 300)
devVzdd_anchors$pmap
dev.off()


png("../../figures/zdVzdd_anchors_homer.png", units = "in", height = 3, width = 4, res = 300)
zdVzdd_anchors$pmap
dev.off()





