# Author: Melyssa Mint0



# load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)
# Read Data in ------------------------------------------------------------

bulk_rna = read_tsv("../../results/DiffExp_RNA/GeneExp_data.tsv")
trap_seq_RNA = read_xlsx("../../sequencing_data/Mellen/mmc1.xlsx", sheet = 1, skip = 1)
trap_seq_pc = read_xls("../../sequencing_data/Mellen/mmc2.xls", sheet = 1)
trap_seq_gc = read_xls("../../sequencing_data/Mellen/mmc2.xls", sheet = 2)
trap_seq_bg = read_xls("../../sequencing_data/Mellen/mmc2.xls", sheet = 3)

# looking at expression  -------------------------------------------------

expr = trap_seq_RNA %>% 
  dplyr::select(symbol, starts_with("TRAPSeq")) %>% 
  dplyr::rename(SYMBOL = symbol, "PC" = "TRAPSeq...7" , "GC" = "TRAPSeq...11" , "BG" = "TRAPSeq...15"  ) %>% 
  pivot_longer(values_to = "FPKM", cols = c("PC", "GC", "BG")) 


expr = expr %>% group_by(SYMBOL) %>% 
  dplyr::slice(which.max(FPKM)) %>% 
  dplyr::mutate(max = name) %>% 
  ungroup() %>% 
  full_join(expr) %>% 
  right_join(bulk_rna)

expr %>% 
  dplyr::filter(!is.na(max)) %>% 
  ggplot(aes(y = gene_sig, fill = max)) +
  geom_bar(position = "fill", stat = "count") 
  

# Looking at differential data --------------------------------------------
bind_rows(list(gc = trap_seq_gc, bg = trap_seq_bg, pc = trap_seq_pc), .id = "id") %>% View()

trap_seq_gc %>% dplyr::mutate(Regulation = "GC") %>% dplyr::select(c("Chromosome", "Start", "End", "Strand", "Gene ID", "Gene Symbol", "Regulation")) %>% 
  full_join(trap_seq_pc %>% dplyr::mutate(Regulation = "PC") %>% dplyr::select(c("Chromosome", "Start", "End", "Strand", "Gene ID", "Gene Symbol", "Regulation"))) %>% 
  full_join(trap_seq_bg %>% dplyr::mutate(Regulation = "BG") %>% dplyr::select(c("Chromosome", "Start", "End", "Strand", "Gene ID", "Gene Symbol", "Regulation"))) %>% 
  dplyr::rename(SYMBOL = `Gene Symbol`) %>% 
  right_join(bulk_rna) %>% 
  group_by(gene_sig) %>% 
  dplyr::count(`Regulation`) %>% 
  dplyr::filter(!is.na(Regulation)) %>% 
  ggplot(aes(x = gene_sig, y =n,  fill = Regulation))+
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label=n),position = position_fill(vjust = 0.5), color = "white", fontface = "bold", size = 8) +
  labs( x = "P60/P7 Gene expression", fill = "TRAP-seq\n Regulation")
  
