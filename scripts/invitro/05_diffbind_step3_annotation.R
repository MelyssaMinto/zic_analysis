# Author: Emilano Sotelo & Melyssa Minto
# This scipt annotates the cutnrun peaks 


# Load libraries ----------------------------------------------------------

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)



# Read in data ------------------------------------------------------------
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
output_dir <- "../../results/invitro/diffbind_cutnrun_zic/"
deseq_obj_contrasts <- readRDS("../../results/invitro/diffbind_cutnrun_zic/deseq_obj_contrasts.Rds")
name_contrast <- "7v3"
res <- readRDS(paste0(output_dir,"res_",name_contrast,"_ashr.Rds"))


# Wrangle data ------------------------------------------------------------
# create peak sets based on differential results
peaks_up <- GenomicRanges::granges(SummarizedExperiment::rowRanges(deseq_obj_contrasts[tidyr::replace_na(res$padj,1) < 0.05 &
                                             tidyr::replace_na(res$log2FoldChange,0) > 0]))
peaks_down <- GenomicRanges::granges(SummarizedExperiment::rowRanges(deseq_obj_contrasts[tidyr::replace_na(res$padj,1) < 0.05 &
                                             tidyr::replace_na(res$log2FoldChange,0) < 0]))
peaks_all <- GenomicRanges::granges(SummarizedExperiment::rowRanges(deseq_obj_contrasts))

elementMetadata(peaks_up) <- res[tidyr::replace_na(res$padj,1) < 0.05 &
                                        tidyr::replace_na(res$log2FoldChange,0) > 0,]
elementMetadata(peaks_down) <- res[tidyr::replace_na(res$padj,1) < 0.05 &
                                        tidyr::replace_na(res$log2FoldChange,0) < 0,]

elementMetadata(peaks_all) <- res

## annotate peaks by differential results
peaks_up_anno <- annotatePeak(peak = peaks_up,TxDb = txdb,
                               level = "gene",
                               assignGenomicAnnotation = T,
                               annoDb = "org.Mm.eg.db")
peaks_down_anno <- annotatePeak(peak = peaks_down,
                               TxDb = txdb,
                               level = "gene",
                               assignGenomicAnnotation = T,
                               annoDb = "org.Mm.eg.db")
peaks_all_anno <- annotatePeak(peak = peaks_all,
                               TxDb = txdb,
                               level = "gene",
                               assignGenomicAnnotation = T,
                               annoDb = "org.Mm.eg.db")


# outputting data ---------------------------------------------------------
saveRDS(peaks_up_anno, file = paste0(output_dir,"peaks_",name_contrast,"_up_annotation.Rds"))
saveRDS(peaks_down_anno, file = paste0(output_dir,"peaks_",name_contrast,"_down_annotation.Rds"))
saveRDS(peaks_all_anno, file = paste0(output_dir,"peaks_",name_contrast,"_all_annotation.Rds"))


readr::write_tsv(as.data.frame(peaks_up_anno@anno), file = paste0(output_dir,"peaks_",name_contrast,"_up_annotation.tsv"))
readr::write_tsv(as.data.frame(peaks_down_anno@anno), file = paste0(output_dir,"peaks_",name_contrast,"_down_annotation.tsv"))
readr::write_tsv(as.data.frame(peaks_all_anno@anno), file = paste0(output_dir,"peaks_",name_contrast,"_all_annotation.tsv"))



