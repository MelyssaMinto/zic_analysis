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

sample_files = c("../../results/DiffExp_ZicChIP/P60vP7_DOWN.bed", "../../results/DiffExp_ZicChIP/P60vP7_NS.bed", "../../results/DiffExp_ZicChIP/P60vP7_UP.bed")




# mapping peaks to nearest gene -------------------------------------------

sample_files <- as.list(sample_files)
names(sample_files) <- c("P7", "N.S.", "P60")

zic_anno = lapply(sample_files, annotatePeak, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion=c(-1000, 1000), verbose=FALSE)

zic_anno_dat = bind_rows(list( P7 = zic_anno$P7@annoStat, N.S. = zic_anno$N.S.@annoStat, P60 = zic_anno$P60@annoStat), .id = "id")


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

fimo = fimo %>% 
  group_by(sequence_name, motif_id) %>% 
  dplyr::mutate(motif_alt_id = 1:n()) %>% 
  ungroup()

# > merging overlap and fimo ----------------------------------------------

zic_peak_data = zic_diff %>% 
  dplyr::mutate(sequence_name = paste0(Chr, ":", Start, "-", End)) %>% 
  dplyr::rename_with( ~ paste0( "zic_", .x)) %>% 
  dplyr::rename(sequence_name = zic_sequence_name, zic_sig = zic_zic_sig) %>% 
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
                sequence_name = ifelse(is.na(sequence_name), "no zic peak", sequence_name),
                motif_name = paste0(motif_id, motif_alt_id)) %>% 
  dplyr::select(sequence_name, starts_with("zic"), ends_with("name"), ends_with("sig"), starts_with("motif"), "matched_sequence") %>% 
  distinct() %>% 
  # counting the number of motifs in peak
  group_by(sequence_name) %>% 
  dplyr::mutate(n_motifs_in_peak = n_distinct(motif_name)) %>% 
  dplyr::mutate(n_dnase_in_peak = n_distinct(dnase_sequence_name)) %>% 
  dplyr::mutate(n_k27ac_in_peak = n_distinct(k27ac_sequence_name)) %>% 
  ungroup() %>%  
  dplyr::mutate(n_motifs_in_peak = ifelse(motif_id %in% "no motif", 0, n_motifs_in_peak),
                n_dnase_in_peak = ifelse(dnase_sequence_name %in% "no overlap" | is.na(dnase_sequence_name), 0, n_dnase_in_peak),
                n_k27ac_in_peak = ifelse(k27ac_sequence_name %in% "no overlap" | is.na(k27ac_sequence_name), 0, n_k27ac_in_peak)) 

#summarizing data to a table
zic_peak_table = zic_peak_data %>% 
  group_by(sequence_name) %>% 
  # dplyr::select(-starts_with(n)) %>% 
  # collapses peaks by loop
  summarise( zic_motifs = glue_collapse(unique(motif_name), sep = ", "),
             k27ac_peaks = glue_collapse(unique(k27ac_sequence_name), sep = ", "),
             dnase_peaks = glue_collapse(unique(dnase_sequence_name), sep = ", ")) %>% 
  ungroup() %>%
  left_join(zic_diff %>% dplyr::mutate(sequence_name = paste0(Chr, ":", Start, "-", End))) %>% 
  dplyr::select(-PeakID)


# Writing output -----------------------------------------------------------
write_tsv(zic_peak_data, "../../results/FinalTables/zic_peak_data.txt")
write_tsv(zic_peak_table, "../../results/FinalTables/zic_peak_table.txt")
write_tsv(zic_anno_dat, "../../results/FinalTables/zic_annotations.txt")

