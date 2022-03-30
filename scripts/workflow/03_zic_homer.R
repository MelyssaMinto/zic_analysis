# Author: Melyssa Minto
# this script will get the trasncriptional enrichment of TFs whose motifs were enriched in 
# early and late Zic peaks

# Load Libraries ----------------------------------------------------------
library(tidyverse)
library(readxl)
library(reshape2)
library(foreach)
library(xlsx)

# Read in data ------------------------------------------------------------

# peak files
peakFilesOfInterest = list.files("../../results/DiffExp_ZicChIP/", pattern = "P60vP7")

# motif enrichment results
homerResultsPath = paste0("../../results//homer_results/", gsub(".bed", "", peakFilesOfInterest) , "/results/knownResults.txt")
homerResults <- homerResultsPath %>%
  map(read_tsv) 


# DiffExp data
DEResults <- read_tsv("../../results/DiffExp_RNA/GeneExp_data.tsv")

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


write_tsv(as.data.frame(homerResultsClean[[1]]),paste0("../../results/homer_results/", names(homerResultsClean)[1], "/homer_results.tsv"))
write_tsv(as.data.frame(homerResultsClean[[2]]),paste0("../../results/homer_results/", names(homerResultsClean)[2], "/homer_results.tsv"))
write_tsv(as.data.frame(homerResultsClean[[3]]),paste0("../../results/homer_results/", names(homerResultsClean)[3], "/homer_results.tsv"))

