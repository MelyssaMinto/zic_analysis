# Author: Emiliano Sotelo & Melyssa Minto
# This scripts formats the cutnrun peaks for differential analysis


# load packages -----------------------------------------------------------
library(rtracklayer)
library(tidyverse)
library(csaw)
library(DESeq2)

# Read in data ------------------------------------------------------------
## load cutnrun peaks 

cutnrun_peaks_dir <- "../../sequencing_data/invitro_data/cutnrun_zic/seacr_peaks/"
cutnrun_peaks_files <- list.files(cutnrun_peaks_dir,pattern = ".bed")
cutnrun_peaks_full <- paste0(cutnrun_peaks_dir,cutnrun_peaks_files)

## alignment files
bam_files <- list.files("../../sequencing_data/invitro_data/cutnrun_zic/bam/",pattern = ".bam$",full.names = T)

mm10_blacklist <- import("../../sequencing_data/invitro_data/bed/mm10-blacklist.v2.bed.gz")


# Wrangle data ------------------------------------------------------------
## cutnrun peaks as granges objects
cutnrun_peaks_granges <- lapply(cutnrun_peaks_full, import.bedGraph)

names(cutnrun_peaks_granges) <- str_split(cutnrun_peaks_files,"_",n=2) %>% map_chr(1)

names(cutnrun_peaks_granges) <- paste0("cutnrun-",names(cutnrun_peaks_granges))

## get union
cutnrun_peaks_union_granges <- GenomicRanges::union(unlist(GRangesList(cutnrun_peaks_granges)),
                                                    unlist(GRangesList(cutnrun_peaks_granges)))
## remove peaks that overlap blacklist
cutnrun_peaks_union_granges <- cutnrun_peaks_union_granges[!overlapsAny(cutnrun_peaks_union_granges, mm10_blacklist)]

## Export union
export.bed(cutnrun_peaks_union_granges,"../../results/invitro/diffbind_cutnrun_zic/cutnrun_zic_seacr_peaks_union.bed")


# Extract counts ----------------------------------------------------------
## extract counts
counts <- csaw::regionCounts(bam.files = bam_files,
                             regions = cutnrun_peaks_union_granges,
                             param = csaw::readParam(dedup = FALSE, minq = 255, pe = "both"))

# format counts matrix
colnames(counts) <- str_split(basename(bam_files),"_",2) %>% map_chr(1)
colData(counts)$rep <- str_split(colnames(counts),"-") %>% map_chr(1)
colData(counts)$time <- str_split(colnames(counts),"-") %>% map_chr(2)


# Set-up differential analysis ---------------------------------------------
# Get DESEq object
deseq_obj <- DESeqDataSet(counts, design = ~ time)
saveRDS(deseq_obj, "../../results/invitro/diffbind_cutnrun_zic/deseq_obj.Rds")

