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
  dplyr::mutate(loop_id = paste0("loop_", 1:n())) 

HiC_Gene = GSE164360_Hi_C_Summary %>% 
  dplyr::select(ends_with(c("A", "B")), GeneID, Ensemble_ID) %>% 
  dplyr::mutate(loop_id = paste0("loop_", 1:n())) 


# reformat data for bedfile
EP_loops = bind_rows(
  list( 
    promoter = EPG_loops %>% dplyr::select(loop_id, starts_with(c("promoter"))) %>% rename_with( ~gsub("promoter_", "", .x), starts_with("promoter")),
    enhancer = EPG_loops %>% dplyr::select(loop_id, starts_with(c("enhancer"))) %>% rename_with( ~gsub("enhancer_", "", .x), starts_with("enhancer"))
  ),
  .id = "id"
) %>% 
  dplyr::arrange(loop_id) %>% 
  # making the anchor bins 10KB 
  dplyr::mutate( mid = end - start/2,
                 new_start = mid - 5000,
                 new_end = mid + 5000)


HiC_Gene_loops = bind_rows(list( A = HiC_Gene %>% 
                  dplyr::select(ends_with("A"), loop_id) %>% 
                  set_names(gsub("A|#", "", names(.))),
                B = HiC_Gene %>% 
                  dplyr::select(ends_with("B"), loop_id) %>% 
                  set_names(gsub("B", "", names(.)))
                  )) %>% 
  dplyr::mutate(chr = ifelse(is.na(chr), "chrX", paste0("chr",chr)))


#  writing output ---------------------------------------------------------


EP_loops %>% 
  dplyr::select(chr, start, end) %>% 
  write_tsv("../../sequencing_data/Reddy/P4_EP_loops.bed")

EP_loops %>% 
  dplyr::select(chr, new_start, new_end) %>% 
  write_tsv("../../sequencing_data/Reddy/P4_EP_loops_ext.bed")

EPG_loops  %>% dplyr::select(loop_id, GeneID, Ensemble_ID) %>% 
  left_join(EP_loops) %>% 
  write_tsv("../../sequencing_data/Reddy/E-P-Gene_map.txt")

HiC_Gene_loops %>% 
  dplyr::select(chr, start, end) %>% 
  write_tsv("../../sequencing_data/Reddy/P4_loops.bed")

HiC_Gene_loops %>% 
  left_join(HiC_Gene %>% dplyr::select(loop_id, GeneID, Ensemble_ID)) %>% 
  write_tsv("../../sequencing_data/Reddy/loops_gene_map.txt")



