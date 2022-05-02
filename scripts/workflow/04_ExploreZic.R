# Author: Melyssa Minto
# This script will explore the zic binding profile in the developing cerebellum


# Load Libraries ----------------------------------------------------------
library(tidyverse)

# Read in data ------------------------------------------------------------
fimo <- read_delim("../../results/FIMO/Zic/fimo.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE, comment = "#")

zic_diff <- read_tsv("../../results/DiffExp_ZicChIP/ZicChIPDA_data.tsv")
dnase_diff <- read_tsv("../../results/DiffExp_DNase/DNaseDA_data.tsv")
k27ac_diff <- read_tsv("../../results/DiffExp_H3K27ac/K27ac_data.tsv")

zic_dnase_overlap <- read_tsv("../../results/mergedPeaks/DNase_zic_overlap.bed", col_names = F)
zic_k27ac_overlap <- read_tsv("../../results/mergedPeaks/H3K27ac_zic_overlap.bed", col_names = F)





# merging data ------------------------------------------------------------


# > Formatting data sets for merging --------------------------------------
zic_dnase_overlap = zic_dnase_overlap %>% 
  dplyr::mutate(sequence_name = paste0(X1, ":", X2, "-", X3),
                dnase_sequence_name = paste0(X4, ":", X5, "-", X6)) %>% 
  dplyr::select(ends_with("name"))

zic_k27ac_overlap = zic_k27ac_overlap %>% 
  dplyr::mutate(sequence_name = paste0(X1, ":", X2, "-", X3),
                k27ac_sequence_name = paste0(X4, ":", X5, "-", X6)) %>% 
  dplyr::select(ends_with("name"))

dnase_diff = dnase_diff %>%  
  dplyr::mutate(dnase_sequence_name = paste0(Chr, ":", Start, "-", End))
  
k27ac_diff = k27ac_diff %>%  
  dplyr::mutate(k27ac_sequence_name = paste0(Chr, ":", Start, "-", End))


# > merging overlap and fimo ----------------------------------------------

zic_peak_data = zic_diff %>% 
  dplyr::mutate(sequence_name = paste0(Chr, ":", Start, "-", End)) %>% 
  dplyr::rename_with( ~ paste0( "zic_", .x)) %>% 
  dplyr::rename(sequence_name = zic_sequence_name, zic_zic_sig = zic_sig) %>% 
  # joining fimo data
  left_join(fimo, by = "sequence_name") %>% 
  # annotating peaks with no motifs
  dplyr::mutate(motif_id = ifelse(is.na(motif_id), "no motif", motif_id)) %>% 
  # joining DNAse and K27ac overlap 
  left_join(zic_dnase_overlap, by = "sequence_name") %>% 
  left_join(zic_k27ac_overlap, by = "sequence_name") %>% 
  # joining Dnase and K27ac diff signal
  full_join(dnase_diff, by = "dnase_sequence_name") %>% 
  full_join(k27ac_diff, by = "k27ac_sequence_name") %>% 
  # annotating Dnase and K27ac peaks with no zic overlap
  dplyr::mutate(dnase_sequence_name = ifelse(is.na(dnase_sequence_name) & !is.na(sequence_name), "no overlap", dnase_sequence_name),
                k27ac_sequence_name = ifelse(is.na(k27ac_sequence_name) & !is.na(sequence_name), "no overlap", k27ac_sequence_name),
                sequence_name = ifelse(is.na(sequence_name), "no zic peak", sequence_name)) %>% 
  dplyr::select(sequence_name, starts_with("zic"), ends_with("name"), ends_with("sig"), "motif_id") %>% 
  distinct() %>% 
  # counting the number of motifs in peak
  group_by(sequence_name) %>% 
  add_count(motif_id, name="n_motifs_in_peak") %>% 
  add_count(dnase_sequence_name, name="n_dnase_in_peak") %>% 
  add_count(k27ac_sequence_name, name="n_k27ac_in_peak") %>% 
  ungroup() %>%  
  dplyr::mutate(n_motifs_in_peak = ifelse(motif_id %in% "no motif", 0, n_motifs_in_peak),
                n_dnase_in_peak = ifelse(dnase_sequence_name %in% "no overlap" | is.na(dnase_sequence_name), 0, n_dnase_in_peak),
                n_k27ac_in_peak = ifelse(k27ac_sequence_name %in% "no overlap" | is.na(k27ac_sequence_name), 0, n_k27ac_in_peak)) 

#summarizing data to a table
zic_peak_table = zic_peak_data %>% 
  group_by(sequence_name) %>% 
  # dplyr::select(-starts_with(n)) %>% 
  # collapses peaks by loop
  summarise( zic_motifs = glue_collapse(unique(motif_id), sep = ", "),
             k27ac_peaks = glue_collapse(unique(k27ac_sequence_name), sep = ", "),
             dnase_peaks = glue_collapse(unique(dnase_sequence_name), sep = ", ")) %>% 
  ungroup() %>%
  left_join(zic_diff %>% dplyr::mutate(sequence_name = paste0(Chr, ":", Start, "-", End))) %>% 
  dplyr::select(-PeakID)


# Writing output -----------------------------------------------------------
write_tsv(zic_peak_data, "../../results/FinalTables/zic_peak_data.txt")
write_tsv(zic_peak_table, "../../results/FinalTables/zic_peak_table.txt")



