# Author: Emilaino Sotelo & Melyssa Minto
# Comparing cut and run peaks to chip peaks


# Load packages -----------------------------------------------------------
library(rtracklayer)
#library(tidyverse)
library(csaw)
library(DESeq2)
library(GenomicRanges)
# load cutnrun peaks as granges objects -----------------------------------
# peaks_files <- list.files("../../sequencing_data/invitro_data/cutnrun_zic/seacr_peaks",pattern = ".bed",full.names = TRUE)
# 
# peaks_files <- c(peaks_files, list.files("../../sequencing_data/invitro_data/cutnrun_zic/chip_zic_seacr_peaks",pattern = ".bed",full.names = TRUE)[1:4])
# 
# peaks_granges <- lapply(peaks_files, import.bedGraph)
# 
# 
# ## get union
# peaks_union_granges <- GenomicRanges::union(unlist(GRangesList(peaks_granges)),
#                                                     unlist(GRangesList(peaks_granges)))
# ## remove peaks that overlap blacklist
# 
# mm10_blacklist <- import("../../../../../genomeData/mm10/bed/mm10-blacklist.bed")
# peaks_union_granges <- peaks_union_granges[!overlapsAny(peaks_union_granges, mm10_blacklist)]

# read in peakset 
peaks_union_granges <- readr::read_tsv("../../sequencing_data/invitro_data/cutnrun_zic/bed/cutnrun_chip_zic_seacr_peaks_union.bed", col_names = F)


# Generate counts ---------------------------------------------------------
# make list of bam files
bam_files <- list.files("../../sequencing_data/invitro_data/cutnrun_zic/bam/",pattern = ".bam$",full.names = T)
bam_files <- c(bam_files, list.files("../../sequencing_data/preprocessed_data/zic_chip/",pattern = ".bam$",full.names = T, recursive = T))


counts_cutnrun <- csaw::regionCounts(bam.files = bam_files[1:12],
                             regions = peaks_union_granges,
                             param = csaw::readParam(dedup = FALSE, minq = 255,max.frag=1000, pe = "both"))

counts_chip <- csaw::regionCounts(bam.files = bam_files[13:16],
                                     regions = peaks_union_granges,
                                     param = csaw::readParam(dedup = FALSE, minq = 255,max.frag=1000))
counts <- cbind(counts_cutnrun, counts_chip)
# Get counts matrix
bam_names <- str_split(basename(bam_files),"_",2) %>% map_chr(1)
bam_names <- str_remove(bam_names,"Aligned.sortedByCoord.out.bam")

bam_names[13:16] <- paste(bam_names[13:16],c("P7","P7","P60","P60"), sep = "-")

colnames(counts) <- bam_names

colData(counts)$rep <- str_split(colnames(counts),"-") %>% map_chr(1)
colData(counts)$time <- str_split(colnames(counts),"-") %>% map_chr(2)

# Create DESeq object -----------------------------------------------------
deseq_obj <- DESeqDataSet(counts, design = ~ time)
saveRDS(deseq_obj, "../../results/invitro/diffbind_cutnrun_zic/deseq_obj_chip_cutnrun.Rds")
