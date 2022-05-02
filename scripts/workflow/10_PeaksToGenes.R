# Author: Melyssa Minto

# Mapping genes to peaks via chromatin loops.


# Load packages -----------------------------------------------------------
library(tidyverse)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(EnsDb.Mmusculus.v79)
library(forcats)
library(glue)
options(scipen=999)

# Read in data ------------------------------------------------------------
# Differential RNA-seq, Zic ChIP-seq, DNase-seq, H3K27ac ChIP-seq results
bulkRNA_DE_data <- read_delim("../../results/DiffExp_RNA/GeneExp_data.tsv",
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)
zic_DA_data <- read_delim("../../results/DiffExp_ZicChIP/ZicChIPDA_data.tsv", 
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)
dnase_DA_data <- read_delim("../../results/DiffExp_DNase/DNaseDA_data.tsv",
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)
k27ac_DA_data <- read_delim("../../results/DiffExp_H3K27ac/K27ac_data.tsv",
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)
# looping data 
# > P56 PLAC-seq loops. source: https://www.nature.com/articles/s41586-019-1190-7
loop_data_adult <- read_delim("../../sequencing_data/Yamada/combined_MAPS_peaks.txt",
                        delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)
# > P4 loops derived from Hi-C. source: https://www.nature.com/articles/s41467-021-25846-3
loop_data_young <- read_delim("../../sequencing_data/Reddy/loops_gene_map.txt",
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)
# > bed file of all anchor regions
anchors <- read_delim("../../results/mergedPeaks/anchors.bed",
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)
# > intersection of peaks and anchors 
zic_loop_map <- read_delim("../../results/mergedPeaks/zic_anchors.bed",
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)
dnase_loop_map <- read_delim("../../results/mergedPeaks/DNase_anchors.bed",
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)
k27ac_loop_map <- read_delim("../../results/mergedPeaks/H3K27ac_anchors.bed",
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)


# Define Functions --------------------------------------------------------
# > function to filter for zic and gene regations and output in bed6 format
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

# Wrangle  data ------------------------------------------------------------
# > wrangling loop data  --------------------------------------------------------
# reshaping loop anchors
loop_data_adult = loop_data_adult %>% dplyr::mutate(loop_id = paste0("loop_", 1:n())) 

loops = bind_rows( 
  list(
    adult = bind_rows(
      loop_data_adult %>% 
        dplyr::select(loop_id, X1:X3) %>% 
        dplyr::rename("chr" = "X1", "start" = "X2", "end"= "X3"),
      loop_data_adult %>% 
        dplyr::select(loop_id, X4:X6) %>% 
        dplyr::rename("chr" = "X4", "start" = "X5", "end"= "X6") ) %>% 
        distinct() ,
  young = loop_data_young %>% 
    dplyr::select(loop_id, chr, start, end)
  ),
  .id = "id")  %>% 
  dplyr::mutate(id = paste0(id,"_", loop_id), 
                anchor_name = paste0(chr, ":", start, "-", end)) %>% 
  dplyr::select(-loop_id) %>% 
  distinct()

  
# > Assigning anchors to genes --------------------------------------------

edb <- EnsDb.Mmusculus.v79
tx <- transcripts(edb, columns=c("tx_id", "gene_id", "gene_name"))
mapping <- data.frame(tx_id=tx$tx_id, SYMBOL=tx$gene_name)
## nearest gene method
annot = as.data.frame(annotatePeak(makeGRangesFromDataFrame(df = loops), TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)@anno) %>% 
  dplyr::mutate(seqstart = start, seqstop = end) %>% 
  dplyr::select(seqnames, seqstart, seqstop, annotation, transcriptId, distanceToTSS) %>% 
  # adding gene symbols
  dplyr::mutate(tx_id = gsub("\\..*","",transcriptId)) %>% 
  left_join(mapping, by = "tx_id" ) %>% 
  dplyr::rename(gene_name = SYMBOL) %>% 
  # creating anchor identifier
  dplyr::mutate(anchor_name = paste0(seqnames, ":", seqstart,'-', seqstop)) %>% 
  full_join(loops, by = "anchor_name") %>% 
  dplyr::select(-seqnames, -seqstart, -seqstop) %>% 
  distinct() %>% 
  # selecting the nearest gene to anchor mapping for each loop
  group_by(id) %>% 
  slice_min(abs(distanceToTSS))  %>% 
  # adding back full list of anchors was not selected as having the closest gene mapping
  full_join(loops %>% dplyr::select(id, anchor_name), by = "id") %>% 
  distinct() %>% 
  dplyr::filter(anchor_name.x != anchor_name.y)
  
head(annot)

## % overlap with gene


# > assign genes to peaks via loops  --------------------------------------
# >>  Zic -----------------------------------------------------------------
zic_data = zic_loop_map %>% 
  dplyr::mutate(anchor_name = paste0(X1, ":", X2, "-", X3),
                zic_peak = paste0(X4, ":", X5, "-", X6)) %>% 
  dplyr::select(anchor_name, zic_peak) %>% 
  # adding Zic regulation
  full_join( zic_DA_data %>% dplyr::mutate(zic_peak = paste0(Chr, ":", Start, "-", End)) %>% dplyr::select(zic_peak, zic_sig)) %>% 
  distinct() %>% 
  # adding loop id
  left_join(loops)
# >> DNase ----------------------------------------------------------------
dnase_data = dnase_loop_map %>% 
  dplyr::mutate(anchor_name = paste0(X1, ":", X2, "-", X3),
                dnase_peak = paste0(X4, ":", X5, "-", X6)) %>% 
  dplyr::select(anchor_name, dnase_peak) %>% 
  # adding DNase regulation
  full_join(dnase_DA_data %>% dplyr::mutate(dnase_peak = paste0(Chr, ":", Start, "-", End)) %>% dplyr::select(dnase_peak, dnase_sig)) %>% 
  distinct() %>% 
  # adding loop id
  left_join(loops)

# >> H3K27ac --------------------------------------------------------------
k27ac_data = k27ac_loop_map %>% 
  dplyr::mutate(anchor_name = paste0(X1, ":", X2, "-", X3),
                k27ac_peak = paste0(X4, ":", X5, "-", X6)) %>% 
  dplyr::select(anchor_name, k27ac_peak) %>% 
  # adding H3K27ac regulation
  full_join(k27ac_DA_data %>% dplyr::mutate(k27ac_peak = paste0(Chr, ":", Start, "-", End)) %>% dplyr::select(k27ac_peak, k27ac_sig)) %>% 
  distinct() %>% 
  # adding loop id
  left_join(loops)


  




# Combining all data -----------------------------------------------------

mapped_data = annot %>% 
  # adding loop maps
  full_join(zic_data, by = "id",na_matches = "never") %>% 
  distinct() %>%
  full_join(dnase_data, by = "id", na_matches = "never") %>% 
  distinct() %>% 
  full_join(k27ac_data,by = "id", na_matches = "never") %>% 
  distinct() %>% 
  # removing anchor identifier
  dplyr::select(-starts_with("anchor"), -tx_id) %>% 
  distinct() %>% 
  # adding bulk RNA data
  full_join(bulkRNA_DE_data %>% dplyr::rename(gene_name = SYMBOL, gene_baseMean = baseMean, gene_padj = padj, gene_lfc = log2FoldChange), by = "gene_name",  na_matches = "never") %>% 
  distinct() %>% 
  #formatting
  dplyr::select(starts_with(c("id","zic", "k27ac", "dnase", "gene", "p7", "p60"))) %>% 
  relocate(id, zic_peak, k27ac_peak, dnase_peak, gene_name) %>% 
  ungroup()

 

mapped_data_table = mapped_data %>% 
  dplyr::select(id, zic_peak, k27ac_peak, dnase_peak) %>% 
  distinct() %>% 
  dplyr::filter(!(id == "no loop")) %>% 
  group_by(id) %>% 
  # collapses peaks by loop
  summarise(zic_peaks = glue_collapse(unique(zic_peak), sep = ", ", ),
            k27ac_peaks = glue_collapse(unique(k27ac_peak), sep = ", "),
            dnase_peaks = glue_collapse(unique(dnase_peak), sep = ", ")) %>% 
  ungroup() %>% 
  # adding loop to gene info
  left_join(mapped_data %>% dplyr::select(id, ends_with("mean"), starts_with("gene") ) %>% distinct(), by = "id") %>% 
  #formatting
  relocate(id, zic_peaks, k27ac_peaks, dnase_peaks, gene_name) 



# Output Mapped Peaks -----------------------------------------------------

write_tsv(mapped_data, "../../results/FinalTables/mapped_data.txt")
write_tsv(mapped_data_table, "../../results/FinalTables/mapped_data_table.txt")

# Output peak sets --------------------------------------------------------
output_peak_set("UP", "UP") %>% 
  write_tsv("../../results/peak_gene/late_activating/P60_peaks_UpGenes.bed", col_names = F)

output_peak_set("UP", "DOWN") %>% 
  write_tsv("../../results/peak_gene/late_repressive/P60_peaks_DOWNGenes.bed", col_names = F)

output_peak_set("DOWN", "DOWN") %>% 
  write_tsv("../../results/peak_gene/early_activating/P7_peaks_DOWNGenes.bed", col_names = F)

output_peak_set("DOWN", "UP") %>% 
  write_tsv("../../results/peak_gene/early_repressive/P7_peaks_UpGenes.bed", col_names = F)

bind_rows(output_peak_set("UP", "UP"),output_peak_set("DOWN", "DOWN") ) %>%
  write_tsv("../../results/peak_gene/activating/activating.bed", col_names = F)

bind_rows(output_peak_set("UP", "DOWN"),output_peak_set("DOWN", "UP") ) %>%
  write_tsv("../../results/peak_gene/repressive/repressive.bed", col_names = F)

output_peak_set("UP", c("UP", "DOWN")) %>% 
  write_tsv("../../results/peak_gene/late/late.bed", col_names = F)

output_peak_set("DOWN", c("UP", "DOWN")) %>% 
  write_tsv("../../results/peak_gene/early/early.bed", col_names = F)

output_peak_set(c("UP", "DOWN", "N.S"), c("UP")) %>% 
  write_tsv("../../results/peak_gene/zicPeaks_P60genes.bed", col_names = F)

output_peak_set(c("UP", "DOWN", "N.S."), c("DOWN")) %>% 
  write_tsv("../../results/peak_gene/zicPeaks_P7genes.bed", col_names = F)

output_peak_set(c("UP", "DOWN", "N.S."), c("N.S.")) %>% 
  write_tsv("../../results/peak_gene/zicPeaks_staticgenes.bed", col_names = F)


