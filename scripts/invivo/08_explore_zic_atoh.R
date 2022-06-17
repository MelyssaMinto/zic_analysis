# Author: Melyssa Minto
# Explore zic, atoh1 overlap


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(EnsDb.Mmusculus.v79)

# Read in data ------------------------------------------------------------
P60vP7_zic = read_tsv("../../results/invivo/DiffExp_ZicChIP/ZicChIPDA_data.tsv")

atoh1_zic_overlap_p7 = read_tsv("../../results/invivo/mergedPeaks/zic_atoh1_p7.bed", col_names = F) 
atoh1_zic_overlap_p60 = read_tsv("../../results/invivo/mergedPeaks/zic_atoh1_p60.bed", col_names = F) 
atoh1_zic_overlap_static = read_tsv("../../results/invivo/mergedPeaks/zic_atoh1_static.bed", col_names = F) 

atoh_peaks =  read_tsv("../../results/invivo/mergedPeaks/atoh1.bed", col_names = F)


# wrangling data ----------------------------------------------------------


atoh1_zic_overlap_p7

# mapping to genes
edb <- EnsDb.Mmusculus.v79
tx <- transcripts(edb, columns=c("tx_id", "gene_id", "gene_name"))
mapping <- data.frame(tx_id=tx$tx_id, SYMBOL=tx$gene_name)
annot = annotatePeak(makeGRangesFromDataFrame(df = P60vP7_zic %>%  dplyr::select(Chr, Start, End)), TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)

# binding all data
P60vP7_zic_overlap = 
  # bind the peaks that overlap at atoh1 at developmentally classified peaks 
  bind_rows(list( p7 = atoh1_zic_overlap_p7 %>% dplyr::mutate(Overlap = T) %>% dplyr::select(Chr= X1,Start =  X2, End =X3, Overlap),
          p60 = atoh1_zic_overlap_p60 %>% dplyr::mutate(Overlap = T) %>% dplyr::select(Chr= X1,Start =  X2, End =X3, Overlap),
          static = atoh1_zic_overlap_static %>% dplyr::mutate(Overlap = T) %>% dplyr::select(Chr= X1,Start =  X2, End =X3, Overlap)), 
          .id = "id") %>% 
  # joining the differential binding stats
  full_join(P60vP7_zic) %>% 
  dplyr::mutate(Overlap = ifelse(is.na(Overlap), FALSE, Overlap)) %>% 
  # joining the annotation stats 
  full_join(as.data.frame(annot@anno) %>% dplyr::rename(Chr = seqnames, Start = start, End = end) ) %>% 
  # pull the annotation data 
  dplyr::mutate(group = str_extract(annotation,"[^\\()]+")) %>% 
  dplyr::select(Chr, Start, End, zic_sig, annotation, group, Overlap)



# Output ------------------------------------------------------------------

P60vP7_zic_overlap %>% write_tsv("../../results/invivo/peak_gene/zic_atoh_overlap.txt")


          
