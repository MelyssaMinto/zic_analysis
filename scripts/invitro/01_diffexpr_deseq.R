# Author: Emiliano Sotelo & Melyssa Minto
# This script performs differential expression analysis for invitro rna-seq data

# load libraries ----------------------------------------------------------
library(tidyverse)
library(DESeq2)
# cite ashr


# Read in data ------------------------------------------------------------
# Read and join counts tables
counts_list <- map(list.files("../../sequencing_data/invitro_data/rnaseq/counts",pattern = "_gene.reads", full.names = TRUE),
                   ~ {
                     counts <- read_tsv(.x, col_names = FALSE)
                     col_name <- basename(.x)
                     col_name <- str_remove(col_name,"_gene.reads")
                     colnames(counts) <- c("gene_names",col_name)
                     return(counts)
                   })


# sample info
sample_info <- read_tsv("../../sequencing_data/invitro_data/rnaseq/rnaseq_kd_filereport_read_run_PRJNA259371_tsv.txt")


# Wrangle data ------------------------------------------------------------
# formatting count data
counts <- purrr::reduce(counts_list, left_join, by = "gene_names") %>% 
  # removing gene counting stats
  dplyr::filter(!grepl("__", gene_names)) %>% 
  # add gene names as rownames
  column_to_rownames("gene_names")



# formatting sample info
sample_info = sample_info %>% 
  dplyr::filter(run_accession %in% colnames(counts)) %>% 
  dplyr::mutate(day = str_split(sample_title, " ",2) %>% map_chr(1),
         repn = str_split(sample_title, " ") %>% map_chr(~ .x[length(.x)]),
         group = case_when(
           day %in% c("d0","d3") ~ "WT",
           str_detect(sample_title,"Zic2 KD") ~ "Zic2-KD" ,
           str_detect(sample_title, "Zic1 KD") ~ "Zic1-KD" ,
           str_detect(sample_title, "KD Cntrl") ~ "Cntrl-KD",
           str_detect(sample_title, "BDNF") ~ "BDNF",
           str_detect(sample_title, "d7 CGN") ~ "WT"
         )) %>%
  dplyr::mutate(group = paste0(day,"_",group),
         repn = paste0("rep",repn),
         run_accession = run_accession,
         col_name = paste(group,repn,sep = "_")) %>% 
  column_to_rownames("col_name")


colnames(counts) <- sample_info$col_name


# Run differential analysis -----------------------------------------------
#> Make DESEQ objects -----------------------------------------------------

dds_all <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = sample_info,
                                  design = ~ group)



dds_zic1 <- dds_all[,str_detect(colnames(dds_all),"d7_Zic1") | str_detect(colnames(dds_all),"d7_Cntrl") ]
dds_zic1$group <- fct_drop(dds_zic1$group)

dds_zic2 <- dds_all[,str_detect(colnames(dds_all),"d7_Zic2") | str_detect(colnames(dds_all),"d7_Cntrl") ]
dds_zic2$group <- fct_drop(dds_zic2$group)

dds_d7_d0 <- dds_all[,str_detect(colnames(dds_all),"d7_WT") | str_detect(colnames(dds_all),"d0_WT") ]
dds_d7_d0$group <- fct_drop(dds_d7_d0$group)

dds_d7_d3 <- dds_all[,str_detect(colnames(dds_all),"d7_WT") | str_detect(colnames(dds_all),"d3_WT")]
dds_d7_d3$group <- fct_drop(dds_d7_d3$group)


# > Test ------------------------------------------------------------------

# Prefiltering
keep <- rowSums(counts(dds_zic1)) >= 10
dds_zic1 <- dds_zic1[keep,]

keep <- rowSums(counts(dds_zic2)) >= 10
dds_zic2 <- dds_zic2[keep,]

keep <- rowSums(counts(dds_d7_d0)) >= 10
dds_d7_d0 <- dds_d7_d0[keep,]

keep <- rowSums(counts(dds_d7_d3)) >= 10
dds_d7_d3 <- dds_d7_d3[keep,]


# > Build Model -----------------------------------------------------------
dds_zic1 <- DESeq(dds_zic1)
dds_zic2 <- DESeq(dds_zic2)
dds_d7_d0 <- DESeq(dds_d7_d0)
dds_d7_d3 <- DESeq(dds_d7_d3)


# > Results ---------------------------------------------------------------
res_zic1 <- results(object = dds_zic1, alpha = 0.05)
res_zic2 <- results(object = dds_zic2, alpha = 0.05)
res_d7_d0 <- results(object = dds_d7_d0, alpha = 0.05)
res_d7_d3 <- results(object = dds_d7_d3, alpha = 0.05)

# Shrinkage
res_zic1_ashr <- lfcShrink(dds_zic1,
                           res  = res_zic1,
                           type="ashr")
res_zic2_ashr <- lfcShrink(dds_zic2,
                           res  = res_zic2,
                           type="ashr")
res_d7_d0_ashr <- lfcShrink(dds_d7_d0,
                            res  = res_d7_d0,
                            type="ashr")
res_d7_d3_ashr <- lfcShrink(dds_d7_d3,
                            res  = res_d7_d3,
                            type="ashr")


# Output results ----------------------------------------------------------

saveRDS(dds_zic1, "../../results/invitro/diffexpr/dds_zic1.Rds")
saveRDS(dds_zic2, "../../results/invitro/diffexpr/dds_zic2.Rds")
saveRDS(dds_d7_d0, "../../results/invitro/diffexpr/dds_d7_d0.Rds")
saveRDS(dds_d7_d3, "../../results/invitro/diffexpr/dds_d7_d3.Rds")


saveRDS(res_zic1_ashr, "../../results/invitro/diffexpr/res_zic1_ashr.Rds")
saveRDS(res_zic2_ashr, "../../results/invitro/diffexpr/res_zic2_ashr.Rds")
saveRDS(res_d7_d0_ashr, "../../results/invitro/diffexpr/res_d7_d0_ashr.Rds")
saveRDS(res_d7_d3_ashr, "../../results/invitro/diffexpr/res_d7_d3_ashr.Rds")


as.data.frame(res_zic1_ashr) %>%
  rownames_to_column(var = "gene_names") %>% 
  readr::write_tsv("../../results/invitro/diffexpr/res_zic1_ashr.tsv")


as.data.frame(res_zic2_ashr) %>%
  rownames_to_column(var = "gene_names") %>% 
  readr::write_tsv( "../../results/invitro/diffexpr/res_zic2_ashr.tsv")

as.data.frame(res_d7_d0_ashr) %>%
  rownames_to_column(var = "gene_names") %>% 
  readr::write_tsv("../../results/invitro/diffexpr/res_d7_d0_ashr.tsv")

as.data.frame(res_d7_d3_ashr) %>%
  rownames_to_column(var = "gene_names") %>% 
  write_tsv("../../results/invitro/diffexpr/res_d7_d3_ashr.tsv")



