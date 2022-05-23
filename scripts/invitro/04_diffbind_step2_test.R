# Author: Emiliano Sotelo & Melyssa Minto
# This script run the differential analysis on cutnrun data

# Load packages -----------------------------------------------------------
library(DESeq2)
library(tidyverse)

# Read in data ------------------------------------------------------------
# Load deseq obj
deseq_obj <- readRDS("../../results/invitro/diffbind_cutnrun_zic/deseq_obj.Rds")

# Run analysis ------------------------------------------------------------


# Inputs
contrasts_list <- list(c("time","D3Z","D1Z"),
                       c("time","D5Z","D3Z"),
                       c("time","D7Z","D5Z"),
                       c("time","D7Z","D3Z"),
                       c("time","D7Z","D1Z"))
res_names <- c("res_3v1","res_5v3","res_7v5","res_7v3","res_7v1")


# Prefiltering
keep <- rowSums(counts(deseq_obj)) >= 20
deseq_obj <- deseq_obj[keep,]

# Build model
deseq_obj_contrasts <- DESeq(deseq_obj)

# Test
res_list <- map(contrasts_list, function(x) results(contrast = x, 
                                                       object = deseq_obj_contrasts,
                                                       alpha = 0.05) )
                    
# Test with LRT analogous to anova 
deseq_obj_lrt <- DESeq(deseq_obj,
                       test="LRT",
                       reduced=~1)

res_all <- results(deseq_obj_lrt, alpha = 0.05)


# Shrinkage
res_ashr_list <- map(res_list, function(x) lfcShrink(dds = deseq_obj_contrasts,
                                       res = x,
                                       type = "ashr"))

res_all_ashr <- lfcShrink(deseq_obj_contrasts,
                          res  = res_all,
                          type="ashr")

# Get diff ranges

get_ranges_up <- function(res, deseq, padj = 0.05,lfc = 0){
  granges(rowRanges(deseq)[tidyr::replace_na(res$padj,1) < padj & tidyr::replace_na(res$log2FoldChange,0) > lfc])
}
get_ranges_down <- function(res, deseq, padj = 0.05,lfc = 0){
  granges(rowRanges(deseq)[tidyr::replace_na(res$padj,1) < padj & tidyr::replace_na(res$log2FoldChange,0) < lfc])
}

res_up_list <- map(res_ashr_list, get_ranges_up, deseq = deseq_obj_contrasts)
res_down_list <- map(res_ashr_list, get_ranges_down, deseq = deseq_obj_contrasts)

res_all_any <- granges(rowRanges(deseq_obj_lrt)[tidyr::replace_na(res_all_ashr$padj,1) < 0.05 ])
names(res_ashr_list) <- res_names
res_7v1_up_lfc2 <- get_ranges_up(res_ashr_list$res_7v1, deseq_obj_contrasts, lfc = 2)
res_7v1_down_lfc2 <-  get_ranges_down(res_ashr_list$res_7v1, deseq_obj_contrasts, lfc = -2)



# Output results ----------------------------------------------------------

out_dir <- "../../results/invitro/diffbind_cutnrun_zic/"

# Write results objects

names(res_ashr_list) <- paste0(out_dir,res_names,"_ashr.Rds")

iwalk(res_ashr_list, ~ saveRDS(.x, .y))
saveRDS(res_all_ashr, paste0(out_dir, "res_all_ashr.Rds"))

# Write results tables
deseq_obj_contrasts_ranges <- as.data.frame(granges(rowRanges(deseq_obj_contrasts)))[,1:3]

names(res_ashr_list) <- paste0(out_dir,res_names,"_ashr.tsv")
iwalk(res_ashr_list, ~ saveRDS(cbind(deseq_obj_contrasts_ranges,as.data.frame(.x)),
                               .y))

readr::write_tsv(cbind(deseq_obj_contrasts_ranges,as.data.frame(res_all_ashr)),
                 paste0(out_dir, "res_all_ashr.tsv"))

# Write deseq fitted objects
saveRDS(deseq_obj_contrasts,  paste0(out_dir, "deseq_obj_contrasts.Rds"))
saveRDS(deseq_obj_lrt,  paste0(out_dir, "deseq_obj_lrt.Rds"))

# DiffBind bedfiles
names(res_up_list) <- paste0(out_dir, res_names,"_up.bed")
iwalk(res_up_list, ~ rtracklayer::export.bed(.x,.y))

names(res_down_list) <- paste0(out_dir, res_names,"_down.bed")
iwalk(res_down_list, ~ rtracklayer::export.bed(.x,.y))

rtracklayer::export.bed(res_all_any, paste0(out_dir, "res_all_any.bed") )
rtracklayer::export.bed(res_7v1_up_lfc2, paste0(out_dir, "res_7v1_up_lfc2.bed") )
rtracklayer::export.bed(res_7v1_down_lfc2, paste0(out_dir, "res_7v1_down_lfc2.bed"))



