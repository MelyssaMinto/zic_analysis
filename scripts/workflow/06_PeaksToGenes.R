# Author: Melyssa Minto

# Mapping genes to peaks


# Load packages -----------------------------------------------------------
library(tidyverse)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(EnsDb.Mmusculus.v79)
library(forcats)
library(glue)
options(scipen=999)

# Read in data ------------------------------------------------------------
scRNA_DE_data <- read_delim("../../results/DiffExp_GCscRNA/scRNA_DE_data.tsv", 
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)
bulkRNA_DE_data <- read_delim("../../results/DiffExp_RNA/GeneExp_data.tsv",
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)
zic_DA_data <- read_delim("../../results/DiffExp_ZicChIP/ZicChIPDA_data.tsv", 
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)
dnase_DA_data <- read_delim("../../results/DiffExp_DNase/DNaseDA_data.tsv",
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)
k27ac_DA_data <- read_delim("../../results/DiffExp_H3K27ac/K27ac_data.tsv",
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)
loop_data_adult <- read_delim("../../sequencing_data/Yamada/combined_MAPS_peaks.txt",
                        delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)
loop_data_young <- read_delim("../../sequencing_data/Reddy/E-P-Gene_map.txt",
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)
zic_loop_map_adult <- read_delim("../../results/mergedPeaks/zic_p56anchors.bed",
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)
dnase_loop_map_adult <- read_delim("../../results/mergedPeaks/DNase_p56anchors.bed",
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)
k27ac_loop_map_adult <- read_delim("../../results/mergedPeaks/H3K27ac_p56anchors.bed",
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)
zic_loop_map_young <- read_delim("../../results/mergedPeaks/zic_p4anchors.bed",
                                delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)
dnase_loop_map_young <- read_delim("../../results/mergedPeaks/DNase_p4anchors.bed",
                                   delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)
k27ac_loop_map_young <- read_delim("../../results/mergedPeaks/H3K27ac_p4anchors.bed",
                                   delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)

# Define Functions --------------------------------------------------------

output_peak_set <-function(zic_reg, gene_reg){
  mapped_data %>% 
    dplyr::select(starts_with(c("zic", "gene"))) %>% 
    distinct() %>% 
    dplyr::filter(!is.na(zic_peak)) %>% 
    dplyr::filter(!is.na(gene_name)) %>% 
    dplyr::filter(zic_sig %in% zic_reg, gene_sig %in% gene_reg) %>% 
    dplyr::select(zic_peak) %>%
    dplyr::mutate(chr = str_extract(zic_peak,"[^\\:]+"),
                  start = str_extract(zic_peak, "(?<=:).+(?=-)"),
                  end = str_extract(zic_peak, "\\b\\w+$")) %>% 
    dplyr::select(chr, start, end) %>% 
    distinct() %>% 
    dplyr::mutate(name = ".", score = 1, strand = ".")
  
}

# Wrangle adult data ------------------------------------------------------------
# > map loops to genes ----------------------------------------------------
# creating a unique indefieir for each loop
loop_data_adult = loop_data_adult %>% 
  dplyr::mutate(loop_id = paste0("loop_", 1:n()))

# reshaping loop anchors
adult_anchors = bind_rows(
  loop_data_adult[,c(1:3, 7)] %>% 
  dplyr::rename("Chr" = "X1", "Start" = "X2", "End"= "X3"),
  loop_data_adult[,4:7] %>% 
  dplyr::rename("Chr" = "X4", "Start" = "X5", "End"= "X6")
) %>% distinct() 

adult_anchors_collapsed = adult_anchors %>% 
  dplyr::mutate(anchor = paste0(Chr, ":", Start, "-", End)) %>%  
  dplyr::select(anchor, loop_id)

# mapping to genes
edb <- EnsDb.Mmusculus.v79
tx <- transcripts(edb, columns=c("tx_id", "gene_id", "gene_name"))
mapping <- data.frame(tx_id=tx$tx_id, SYMBOL=tx$gene_name)

annot = as.data.frame(annotatePeak(makeGRangesFromDataFrame(df = adult_anchors), TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)@anno) %>% 
  dplyr::mutate(seqstart = start, seqstop = end) %>% 
  dplyr::select(seqnames, seqstart, seqstop, annotation, transcriptId, distanceToTSS) %>% 
  # adding gene symbols
  dplyr::mutate(tx_id = gsub("\\..*","",transcriptId)) %>% 
  left_join(mapping, by = "tx_id" ) %>% 
  dplyr::rename(gene_name = SYMBOL) %>% 
  # creating anchor identifier
  dplyr::mutate(anchor = paste0(seqnames, ":", seqstart,'-', seqstop)) %>% 
  # selecting relavant data
  dplyr::select(anchor, gene_name,tx_id , annotation, distanceToTSS ) %>% 
  # adding back full list of anchors that were not mapped to gene
  full_join(adult_anchors %>% dplyr::mutate(anchor = paste0(Chr, ":", Start, "-", End)) %>%  dplyr::select(anchor, loop_id) ) %>% 
  distinct() %>% 
  # selecting the nearest gene to anchor mapping for each loop
  group_by(loop_id) %>% 
  slice_min(abs(distanceToTSS)) %>% 
  # adding back full list of anchors was not selected as having the closest gene mapping
  full_join(adult_anchors_collapsed, by = "loop_id") %>% 
  distinct() %>% 
  dplyr::filter(anchor.x != anchor.y)


rm(edb, tx, mapping)

# > map loops to Zic ----------------------------------------------------
zic_data_adult = zic_loop_map_adult %>% 
  dplyr::mutate(anchor = paste0(X1, ":", X2, "-", X3),
                zic_peak = paste0(X4, ":", X5, "-", X6)) %>% 
  dplyr::select(anchor, zic_peak) %>% 
  # adding Zic regulation
  full_join( zic_DA_data %>% dplyr::mutate(zic_peak = paste0(Chr, ":", Start, "-", End)) %>% dplyr::select(zic_peak, zic_sig)) %>% 
  distinct() %>% 
  # adding loop id
  left_join(adult_anchors_collapsed)

# > map loops to DNase --------------------------------------------------
dnase_data_adult = dnase_loop_map_adult %>% 
  dplyr::mutate(anchor = paste0(X1, ":", X2, "-", X3),
                dnase_peak = paste0(X4, ":", X5, "-", X6)) %>% 
  dplyr::select(anchor, dnase_peak) %>% 
  # adding DNase regulation
  full_join(dnase_DA_data %>% dplyr::mutate(dnase_peak = paste0(Chr, ":", Start, "-", End)) %>% dplyr::select(dnase_peak, dnase_sig)) %>% 
  distinct() %>% 
  # adding loop id
  left_join(adult_anchors_collapsed)
# > map loops to  H3K27ac ----------------------------------------------
k27ac_data_adult = k27ac_loop_map_adult %>% 
  dplyr::mutate(anchor = paste0(X1, ":", X2, "-", X3),
                k27ac_peak = paste0(X4, ":", X5, "-", X6)) %>% 
  dplyr::select(anchor, k27ac_peak) %>% 
  # adding H3K27ac regulation
  full_join(k27ac_DA_data %>% dplyr::mutate(k27ac_peak = paste0(Chr, ":", Start, "-", End)) %>% dplyr::select(k27ac_peak, k27ac_sig)) %>% 
  distinct() %>% 
  # adding loop id
  left_join(adult_anchors_collapsed)
  





# Wrangle young data ------------------------------------------------------

EP_loops = bind_rows(
  list( 
    promoter = loop_data_young %>% dplyr::select(loop_id, starts_with(c("promoter"))) %>% rename_with( ~gsub("promoter_", "", .x), starts_with("promoter")),
    enhancer = loop_data_young %>% dplyr::select(loop_id, starts_with(c("enhancer"))) %>% rename_with( ~gsub("enhancer_", "", .x), starts_with("enhancer"))
  ),
  .id = "id"
) %>% 
  dplyr::arrange(loop_id) %>% 
  dplyr::left_join( loop_data_young %>% dplyr::select(loop_id, gene_name=GeneID) %>% distinct()) %>% 
  dplyr::mutate(anchor = paste0(chr, ":", start, "-", end))

# loops are already mapped to genes from the Reddy et. al 2021

# > map loops to Zic ------------------------------------------------------

zic_data_young = zic_loop_map_young %>% 
  dplyr::mutate(anchor = paste0(X1, ":", X2, "-", X3),
                zic_peak = paste0(X4, ":", X5, "-", X6)) %>% 
  dplyr::select(anchor, zic_peak) %>% 
  full_join( zic_DA_data %>% dplyr::mutate(zic_peak = paste0(Chr, ":", Start, "-", End)) %>% dplyr::select(zic_peak, zic_sig)) %>% 
  distinct() %>% 
  # adding loop id
  left_join(EP_loops %>% dplyr::select(loop_id, anchor))
# > map loops to Dnase ----------------------------------------------------
dnase_data_young = dnase_loop_map_young %>% 
  dplyr::mutate(anchor = paste0(X1, ":", X2, "-", X3),
                dnase_peak = paste0(X4, ":", X5, "-", X6)) %>% 
  dplyr::select(anchor, dnase_peak) %>% 
  full_join( dnase_DA_data %>% dplyr::mutate(dnase_peak = paste0(Chr, ":", Start, "-", End)) %>% dplyr::select(dnase_peak, dnase_sig)) %>% 
  distinct() %>% 
  # adding loop id
  left_join(EP_loops %>% dplyr::select(loop_id, anchor))

# > map loops to H3K27ac --------------------------------------------------
k27ac_data_young = k27ac_loop_map_young %>% 
  dplyr::mutate(anchor = paste0(X1, ":", X2, "-", X3),
                k27ac_peak = paste0(X4, ":", X5, "-", X6)) %>% 
  dplyr::select(anchor, k27ac_peak) %>% 
  full_join( k27ac_DA_data %>% dplyr::mutate(k27ac_peak = paste0(Chr, ":", Start, "-", End)) %>% dplyr::select(k27ac_peak, k27ac_sig)) %>% 
  distinct() %>% 
  # adding loop id
  left_join(EP_loops %>% dplyr::select(loop_id, anchor))

# combining all data -----------------------------------------------------
mapped_data_adult = annot %>% 
  # adding loop maps
  full_join(zic_data_adult, by = "loop_id",na_matches = "never") %>% 
  distinct() %>%
  full_join(dnase_data_adult, by = "loop_id", na_matches = "never") %>% 
  distinct() %>% 
  full_join(k27ac_data_adult,by = "loop_id", na_matches = "never") %>% 
  distinct() %>% 
  # removing anchor identifier
  dplyr::select(-starts_with("anchor"), -tx_id) %>% 
  distinct() %>% 
  # adding bulk RNA data
  full_join(bulkRNA_DE_data %>% dplyr::rename(gene_name = SYMBOL, gene_baseMean = baseMean, gene_padj = padj, gene_lfc = log2FoldChange), by = "gene_name",  na_matches = "never") %>% 
  distinct() %>% 
  #formatting
  dplyr::select(starts_with(c("loop","zic", "k27ac", "dnase", "gene", "p7", "p60"))) %>% 
  relocate(loop_id, zic_peak, k27ac_peak, dnase_peak, gene_name) %>% 
  ungroup()


mapped_data_young = EP_loops %>% 
  # adding loop maps
  full_join(zic_data_young, by = "loop_id",na_matches = "never") %>% 
  distinct() %>%
  full_join(dnase_data_young, by = "loop_id", na_matches = "never") %>% 
  distinct() %>% 
  full_join(k27ac_data_young,by = "loop_id", na_matches = "never") %>% 
  distinct() %>% 
  # removing anchor identifier
  dplyr::select(-starts_with("anchor")) %>% 
  distinct() %>% 
  # adding bulk RNA data
  full_join(bulkRNA_DE_data %>% dplyr::rename(gene_name = SYMBOL, gene_baseMean = baseMean, gene_padj = padj, gene_lfc = log2FoldChange), by = "gene_name",  na_matches = "never") %>% 
  distinct() %>% 
  #formatting
  dplyr::select(starts_with(c("loop","zic", "k27ac", "dnase", "gene", "p7", "p60"))) %>% 
  relocate(loop_id, zic_peak, k27ac_peak, dnase_peak, gene_name) %>% 
  ungroup()
  

mapped_data = bind_rows( list ( adult = mapped_data_adult, young = mapped_data_young), .id = "id")

mapped_data_table_adult =  annot %>% 
  # adding loop maps
  left_join(zic_data_adult, by = "loop_id",na_matches = "never") %>% 
  distinct() %>%
  left_join(dnase_data_adult, by = "loop_id", na_matches = "never") %>% 
  distinct() %>% 
  left_join(k27ac_data_adult,by = "loop_id", na_matches = "never") %>% 
  distinct() %>% 
  # selecting peak identifier columns
  dplyr::select(loop_id, zic_peak, k27ac_peak, dnase_peak) %>% 
  distinct() %>% 
  group_by(loop_id) %>% 
  # collapses peaks by loop
  summarise(zic_peaks = glue_collapse(unique(zic_peak), sep = ", ", ),
            k27ac_peaks = glue_collapse(unique(k27ac_peak), sep = ", "),
            dnase_peaks = glue_collapse(unique(dnase_peak), sep = ", ")) %>% 
  ungroup() %>% 
  # adding loop to gene info
  left_join(annot, by = "loop_id") %>% 
  dplyr::select("loop_id", ends_with("peaks"), "gene_name") %>% 
  # adding bulk RNA data
  left_join(bulkRNA_DE_data %>% dplyr::rename(gene_name = SYMBOL, gene_baseMean = baseMean, gene_padj = padj, gene_lfc = log2FoldChange), by = "gene_name",  na_matches = "never") %>% 
  distinct() %>%  
  #formatting
  dplyr::select(starts_with(c("loop","zic", "k27ac", "dnase", "gene", "p7","60"))) %>% 
  relocate(loop_id, zic_peaks, k27ac_peaks, dnase_peaks, gene_name) 

mapped_data_table_young = EP_loops %>% 
  # adding loop maps
  full_join(zic_data_young, by = "loop_id",na_matches = "never") %>% 
  distinct() %>%
  full_join(dnase_data_young, by = "loop_id", na_matches = "never") %>% 
  distinct() %>% 
  full_join(k27ac_data_young,by = "loop_id", na_matches = "never") %>% 
  distinct() %>% 
  # selecting peak identifier columns
  dplyr::select(loop_id, zic_peak, k27ac_peak, dnase_peak) %>% 
  distinct() %>% 
  group_by(loop_id) %>% 
  # collapses peaks by loop
  summarise(zic_peaks = glue_collapse(unique(zic_peak), sep = ", ", ),
            k27ac_peaks = glue_collapse(unique(k27ac_peak), sep = ", "),
            dnase_peaks = glue_collapse(unique(dnase_peak), sep = ", ")) %>% 
  ungroup() %>% 
  # adding loop to gene info
  left_join(annot, by = "loop_id") %>% 
  dplyr::select("loop_id", ends_with("peaks"), "gene_name") %>% 
  # adding bulk RNA data
  left_join(bulkRNA_DE_data %>% dplyr::rename(gene_name = SYMBOL, gene_baseMean = baseMean, gene_padj = padj, gene_lfc = log2FoldChange), by = "gene_name",  na_matches = "never") %>% 
  distinct() %>%  
  #formatting
  dplyr::select(starts_with(c("loop","zic", "k27ac", "dnase", "gene", "p7","60"))) %>% 
  relocate(loop_id, zic_peaks, k27ac_peaks, dnase_peaks, gene_name) 

mapped_data_table = bind_rows( list ( adult = mapped_data_table_adult, young = mapped_data_table_young), .id = "id") %>% 
  dplyr::mutate(loop_id = paste0(id,"_" , loop_id))






# Output Mapped Peaks -----------------------------------------------------

write_tsv(mapped_data, "../../results/FinalTables/mapped_data.txt")
write_tsv(mapped_data_table, "../../results/FinalTables/mapped_data_table.txt")

# Output peak sets --------------------------------------------------------
output_peak_set("UP", "UP") %>% write_tsv("../../results/peak_gene/late_activating/P60_peaks_UpGenes.bed", col_names = F)
output_peak_set("UP", "DOWN") %>% write_tsv("../../results/peak_gene/late_repressive/P60_peaks_DOWNGenes.bed", col_names = F)
output_peak_set("DOWN", "DOWN") %>% write_tsv("../../results/peak_gene/early_activating/P7_peaks_DOWNGenes.bed", col_names = F)
output_peak_set("DOWN", "UP") %>% write_tsv("../../results/peak_gene/early_repressive/P7_peaks_UpGenes.bed", col_names = F)

bind_rows(output_peak_set("UP", "UP"),output_peak_set("DOWN", "DOWN") ) %>% 
  write_tsv("../../results/peak_gene/activating/activating.bed", col_names = F)

bind_rows(output_peak_set("UP", "DOWN"),output_peak_set("DOWN", "UP") ) %>% 
  write_tsv("../../results/peak_gene/repressive/repressive.bed", col_names = F)

output_peak_set("UP", c("UP", "DOWN")) %>% write_tsv("../../results/peak_gene/late/late.bed", col_names = F)

output_peak_set("DOWN", c("UP", "DOWN")) %>% write_tsv("../../results/peak_gene/early/early.bed", col_names = F)

# mapped_data %>% 
#   dplyr::filter(!is.na(gene_name) & is.na(zic_peak) & !is.na(loop_id)) %>% 
#   dplyr::select(gene_name) %>% 
#   distinct()

