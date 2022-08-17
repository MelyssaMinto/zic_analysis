# Author: Emilaino Sotelo & Melyssa Minto
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
# looping data 
# > P56 PLAC-seq loops. source: https://www.nature.com/articles/s41586-019-1190-7
loop_data_adult <- read_delim("../../sequencing_data/Yamada/combined_MAPS_peaks.txt",
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)
# > P4 loops derived from Hi-C. source: https://www.nature.com/articles/s41467-021-25846-3
loop_data_young <- read_delim("../../sequencing_data/Reddy/loops_gene_map.txt",
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)
# > bed file of all anchor regions
anchors <- read_delim("../../results/invitro/diffbind_cutnrun_zic/anchors.bed",
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)
# > intersection of peaks and anchors 
zic_loop_map <- read_delim("../../results/invitro/diffbind_cutnrun_zic/zic_anchors.bed",
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = F)



# Wrangle  data ------------------------------------------------------------
# > wrangling loop data  --------------------------------------------------------
# reshaping loop anchors:
# Adding a loop id to every loop, then separate loop anchor pairs to get an object
# with only one of the anchors
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
  dplyr::filter(anchor_name.x != anchor_name.y) # anchor_name.y is the one mapped to genes
  


# > assign genes to peaks via loops  --------------------------------------
# >>  Zic -----------------------------------------------------------------
zic_data = zic_loop_map %>% 
  dplyr::mutate(anchor_name = paste0(X1, ":", X2, "-", X3),
                zic_peak = paste0(X4, ":", X5, "-", X6)) %>% 
  dplyr::select(anchor_name, zic_peak) %>% 
  # adding loop id
  left_join(loops)


# Combining all data -----------------------------------------------------

mapped_data = annot %>% 
  # adding loop maps
  full_join(zic_data, by = "id",na_matches = "never") %>% 
  distinct() %>%
  # removing anchor identifier
  dplyr::select(-starts_with("anchor"), -tx_id) %>% 
  distinct() %>% 
  #formatting
  dplyr::select(starts_with(c("id","zic", "gene"))) %>% 
  relocate(id, zic_peak, gene_name) %>% 
  ungroup()

mapped_data_with_anchors = annot %>% 
  # adding loop maps
  full_join(zic_data, by = "id",na_matches = "never") %>% 
  relocate(id, zic_peak, gene_name, anchor_name.y, distanceToTSS) %>% 
  distinct()

 
mapped_data_table = mapped_data %>% 
  dplyr::select(id, zic_peak) %>% 
  distinct() %>% 
  dplyr::filter(!(id == "no loop")) %>% 
  group_by(id) %>% 
  # collapses peaks by loop
  summarise(zic_peaks = glue_collapse(unique(zic_peak), sep = ", ", )) %>% 
  ungroup() %>% 
  # adding loop to gene info
  left_join(mapped_data %>% dplyr::select(id, ends_with("mean"), starts_with("gene") ) %>% distinct(), by = "id") %>% 
  #formatting
  relocate(id, zic_peaks, gene_name) 


# Output Mapped Peaks -----------------------------------------------------

write_tsv(mapped_data, "../../results/FinalTables/in_vitro_mapped_data.txt")
write_tsv(mapped_data_table, "../../results/FinalTables/in_vitro_mapped_data_table.txt")


