# Melyssa Minto 
# use Reddy et al P4 Hi-C enhancer and promoter loops 


# Load libraries ----------------------------------------------------------
library(tidyverse)
library(readxl)

# Read in data ------------------------------------------------------------


GSE164360_Hi_C_Summary <- read_excel("/media/west-lab-share/Melyssa_Minto/Cerebellum/zic_analysis/sequencing_data/Reddy/GSE164360_Hi-C_Summary.xlsx")


# Wrangle data ------------------------------------------------------------

# Create a bed file of P4 enhancer promoter anchors

# subset for enhancer, promoter and genes for gene to loop paing
EPG_loops = GSE164360_Hi_C_Summary %>% 
  dplyr::select(ends_with(c("chr", "start", "end")), GeneID, Ensemble_ID) %>% 
  dplyr::mutate(loop_id = paste0("loop_", 1:n())) %>% 
  write_tsv("../../sequencing_data/Reddy/E-P-Gene_map.txt")

# reformat data for bedfile
EP_loops = bind_rows(
  list( 
    promoter = EP_loops %>% dplyr::select(loop_id, starts_with(c("promoter"))) %>% rename_with( ~gsub("promoter_", "", .x), starts_with("promoter")),
    enhancer = EP_loops %>% dplyr::select(loop_id, starts_with(c("enhancer"))) %>% rename_with( ~gsub("enhancer_", "", .x), starts_with("enhancer"))
  ),
  .id = "id"
) %>% 
  dplyr::arrange(loop_id)

EP_loops %>% 
  dplyr::select(chr, start, end) %>% 
  write_tsv("../../sequencing_data/Reddy/P4_EP_loops.bed")



