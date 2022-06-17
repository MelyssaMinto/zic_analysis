
# Author: Melyssa Minto (msm110@duke.edu)
# this script compute the differential signal


# Load libs ---------------------------------------------------------------

library(Rsubread)
library(tidyverse)
library(DESeq2)
library(readxl)
library(foreach)
options("scipen"=0, "digits"=7) 

# importing the sample sheet
SAMPLESHEET <- read_excel("../../sequencing_data/SampleSheet_Frank2015.xlsx")

# Zic  --------------------------------------------------------------------

# > Getting peak reads ----------------------------------------------------
# subset for the Zic data in the cerebellum
sample_sheet = SAMPLESHEET %>% 
  dplyr::filter(Assay_Type  %in% "ChIP-Seq") %>% 
  dplyr::filter(antibody %in% "anti-Zic 1/2 C-terminus") %>% 
  dplyr::mutate(group = ifelse(developmental_stage %in% "Postnatal day 7 (P7)", "P7", "P60"))

# read in merged peak file to get reads from
zic_merged <- read_delim("../../results/mergedPeaks/zic.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
chrom_sizes <- read_delim("../../../../../genomeData/mm10/mm10.chrom.sizes.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE) %>% 
  pivot_wider(names_from = X1, values_from = X2)

# formatting as gtf (geneID, Chr, Start, End, Strand)
peakID = paste0("peak", 1:nrow(zic_merged))

all_Zic_sites = data.frame(GeneID = peakID,
                           Chr=zic_merged$X1,
                           Start=zic_merged$X2,
                           End=zic_merged$X3,
                           Strand = rep("",nrow(zic_merged)))

files = files = list.files("../../sequencing_data/preprocessed_data/zic_chip/", pattern = ".bam$", recursive = T)

Zic_counts = featureCounts(paste0("../../sequencing_data/preprocessed_data/zic_chip/", files), 
                           annot.ext = all_Zic_sites,
                           nthreads = 8)



# > Differential binding analysis -----------------------------------------
# extract count matrix
counts <- as.data.frame(Zic_counts$counts) 
colnames(counts) = c("Zic_P7_1", "Zic_P7_2", "Zc_P60_1", "Zic_P0_2")

#run DESeq2
dds<-DESeqDataSetFromMatrix(countData = counts, 
                            colData = sample_sheet,
                            design = ~ group)
dds<-DESeq(dds)

# format results
P60vP7_zic = results(dds, contrast = c("group", "P60", "P7"), tidy = TRUE) %>% 
  dplyr::rename(GeneID = row) %>% 
  left_join(all_Zic_sites) %>% 
  dplyr::select(-Strand) %>% 
  dplyr::mutate(zic_sig = case_when(log2FoldChange >= 0 & padj < 0.05 ~ "UP",
                                     log2FoldChange <= 0 & padj < 0.05 ~ "DOWN",
                                     is.na(padj) ~ "filtered",
                                     TRUE ~ "N.S.")) %>% 
  dplyr::rename(PeakID = GeneID)

# save results
#> differential expression table
as.data.frame(counts(dds, normalized = T)) %>%
  dplyr::mutate(PeakID = P60vP7_zic$PeakID) %>% 
  dplyr::mutate(p7_mean = rowMeans(dplyr::select(., contains("P7"))),
         p60_mean = rowMeans(dplyr::select(.,contains("P60")))) %>% 
  dplyr::select(PeakID, p7_mean, p60_mean) %>% 
  full_join(P60vP7_zic) %>% write_tsv("../../results/invivo/DiffExp_ZicChIP/ZicChIPDA_data.tsv")


#> bed files for P60 peaks
P60vP7_zic %>% 
  dplyr::filter(zic_sig %in% "UP") %>% 
  dplyr::select(Chr, Start, End) %>% 
  write_tsv("../../results/invivo/DiffExp_ZicChIP/P60vP7_UP.bed",col_names = FALSE)
P60vP7_zic %>% 
  dplyr::filter(zic_sig %in% "UP") %>% 
  dplyr::select(Chr, Start, End) %>% 
  dplyr::mutate(name = ".", score = 1, strand = ".") %>% 
  write_tsv("../../results/invivo/DiffExp_ZicChIP/late/late.bed",col_names = FALSE)

#> bed files for P7 peaks
P60vP7_zic %>% 
  dplyr::filter(zic_sig %in% "DOWN") %>% 
  dplyr::select(Chr, Start, End) %>% 
  write_tsv("../../results/invivo/DiffExp_ZicChIP/P60vP7_DOWN.bed",col_names = FALSE)

P60vP7_zic %>% 
  dplyr::filter(zic_sig %in% "DOWN") %>% 
  dplyr::select(Chr, Start, End) %>% 
  dplyr::mutate(name = ".", score = 1, strand = ".") %>% 
  write_tsv("../../results/invivo/DiffExp_ZicChIP/early/early.bed",col_names = FALSE)


#> bed files for static peaks
P60vP7_zic %>% 
  dplyr::filter(zic_sig %in% "N.S.") %>% 
  dplyr::select(Chr, Start, End) %>% 
  write_tsv("../../results/invivo/DiffExp_ZicChIP/P60vP7_NS.bed",col_names = FALSE)

P60vP7_zic %>% 
  dplyr::filter(zic_sig %in% "N.S.") %>% 
  dplyr::select(Chr, Start, End) %>% 
  dplyr::mutate(name = ".", score = 1, strand = ".") %>% 
  write_tsv("../../invivo/results/DiffExp_ZicChIP/static/static.bed",col_names = FALSE)

# > Plotting --------------------------------------------------------------
nums = table(P60vP7_zic$zic_sig)

P60vP7_zic %>% 
  dplyr::filter(zic_sig != "filtered") %>% 
  ggplot(aes(x = log10(baseMean), y = log2FoldChange, color = zic_sig))+
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values= c("blue", "black", "red")) +
  geom_text(x = 3.5, y = 5, label = paste0("Increased Zic Binding\n n = ",nums["UP"] ), color = "red")+
  geom_text(x = 3.5, y = -5, label = paste0("Decreased Zic Binding\n n = ",nums["DOWN"] ), color = "blue")+
  theme_classic() + 
  labs( x = "log(Avg Signal)", y = "Log2FC")+
  theme(legend.position = "none")

ggsave("../../figures/zic_ma.png")  

P60vP7_zic %>% 
  dplyr::filter(zic_sig != "filtered") %>% 
  ggplot(aes(x = log10(baseMean), y = log2FoldChange, color = zic_sig))+
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values= c("blue", "black", "red")) +
  geom_text(x = 3.5, y = 5, label = paste0("Increased Zic Binding\n n = ",nums["UP"] ), color = "red")+
  geom_text(x = 3.5, y = -5, label = paste0("Decreased Zic Binding\n n = ",nums["DOWN"] ), color = "blue")+
  theme_classic() + 
  xlim(c(.7, 4.7))+
  labs( x = "log10(Avg Zic ChIP  Signal)", y = "Log2FC")+
  theme(legend.position = "none")

ggsave("../../figures/zic_ma_crop.png") 


# DNase -------------------------------------------------------------------

# > Getting the DNase Peak reads ------------------------------------------

# subset for the DHS data in the cerebellum
sample_sheet = SAMPLESHEET %>% 
  dplyr::filter(Assay_Type %in% "DNase-Hypersensitivity") %>% 
  dplyr::filter(source_name %in% "whole cerebellum") %>% 
  dplyr::filter(developmental_stage %in% c("Postnatal day 7 (P7)", "Postnatal day 60 (P60)")) %>% 
  dplyr::mutate(group = ifelse(developmental_stage %in% "Postnatal day 7 (P7)", "P7", "P60"))

# read in merged peak file to get reads from
merged <- read_delim("../../results/invivo/mergedPeaks/DNase.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)


# formatting as gtf (geneID, Chr, Start, End, Strand)
peakID = paste0("peak", 1:nrow(merged))
all_bg_sites = data.frame(GeneID = peakID,
                          Chr=merged$X1,
                          Start=merged$X2,
                          End=merged$X3,
                          Strand = rep("",nrow(merged)))

files = list.files("../../sequencing_data/preprocessed_data/DNase/", pattern = "bowtie.bam", recursive = T)
DNAse_counts = featureCounts(paste0("../../sequencing_data/preprocessed_data/DNase/",files), 
                             annot.ext = all_bg_sites,
                             nthreads = 8)

# > Differential Accessibility Analysis -----------------------------------

# extract count matrix
counts <- as.data.frame(DNAse_counts$counts) 
colnames(counts) = c("DNase_P7_1", "DNase_P7_2", "DNase_P7_3", "DNase_P60_1", "DNase_P60_2", "DNase_6P0_3")
# removing p60 rep 2
counts = counts[,-5]
sample_sheet = sample_sheet[-5,]

# run DESeq
DNase_dds<-DESeqDataSetFromMatrix(countData = counts, 
                                  colData = sample_sheet,
                                  design = ~ group)
DNase_dds<-DESeq(DNase_dds)

# format resutls
P60vP7_DNase = results(DNase_dds, contrast = c("group", "P60", "P7"), tidy = TRUE) %>% 
  dplyr::rename(GeneID = row) %>% 
  left_join(all_bg_sites) %>%  
  dplyr::select(-Strand) %>% 
  dplyr::mutate(dnase_sig = case_when(log2FoldChange >= 0 & padj < 0.05 ~ "UP",
                                                               log2FoldChange <= 0 & padj < 0.05 ~ "DOWN",
                                                               is.na(padj) ~ "filtered",
                                                               TRUE ~ "N.S."))



# save results
P60vP7_DNase %>% write_tsv("../../results/invivo/DiffExp_DNase/DNaseDA_data.tsv")

# > Plotting --------------------------------------------------------------
nums = table(P60vP7_DNase$dnase_sig)

P60vP7_DNase %>% 
  dplyr::filter(dnase_sig != "filtered") %>% 
  ggplot(aes(x = log10(baseMean), y = log2FoldChange, color = dnase_sig))+
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values= c("blue", "black", "red")) +
  geom_text(x = 3, y = 4, label = paste0("Opening DHS \n n = ",nums["UP"] ), color = "red")+
  geom_text(x = 3, y = -2, label = paste0("Closing DHS \n n = ",nums["DOWN"] ), color = "blue")+
  theme_classic() + 
  labs( x = "log(Avg Signal)", y = "Log2FC")+
  theme(legend.position = "none")

ggsave("../../figures/DNase_ma.png")  


P60vP7_DNase %>% 
  dplyr::filter(dnase_sig != "filtered") %>% 
  ggplot(aes(x = log10(baseMean), y = log2FoldChange, color = dnase_sig))+
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values= c("blue", "black", "red")) +
  geom_text(x = 3, y = 4, label = paste0("Opening DHS \n n = ",nums["UP"] ), color = "red")+
  geom_text(x = 3, y = -2, label = paste0("Closing DHS \n n = ",nums["DOWN"] ), color = "blue")+
  theme_classic() + 
  labs( x = "log(Avg Signal)", y = "Log2FC")+
  theme(legend.position = "none") +
  xlim(c(.7, 3.5))

ggsave("../../figures/DNase_ma_crop.png") 



# H3K27ac -----------------------------------------------------------------

# subset for the H3K27ac data in the cerebellum
sample_sheet = SAMPLESHEET %>% 
  dplyr::filter(Assay_Type %in% "ChIP-Seq") %>% 
  dplyr::filter(antibody %in% "H3K27Ac, Abcam Ab4729") %>% 
  dplyr::filter(developmental_stage %in% c("Postnatal day 7 (P7)", "Postnatal day 60 (P60)"))%>% 
  dplyr::mutate(group = ifelse(developmental_stage %in% "Postnatal day 7 (P7)", "P7", "P60"))

# read in merged peak file to get reads from
merged <- read_delim("../../results/invivo/mergedPeaks/H3K27ac.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)


# formatting as gtf (geneID, Chr, Start, End, Strand)
peakID = paste0("peak", 1:nrow(merged))
all_bg_sites = data.frame(GeneID = peakID,
                          Chr=merged$X1,
                          Start=merged$X2,
                          End=merged$X3,
                          Strand = rep("",nrow(merged)))

files = list.files("../../sequencing_data/preprocessed_data/H3K27ac//", pattern = "hits.bam", recursive = T)
K27ac_counts = featureCounts(paste0("../../sequencing_data/preprocessed_data/H3K27ac/",files), 
                             annot.ext = all_bg_sites,
                             nthreads = 8)


# > Getting H3K27ac peak reads --------------------------------------------

# extract count matrix
counts <- as.data.frame(K27ac_counts$counts) 
colnames(counts) = c("H3K27ac_P7_1", "H3K27ac_P7_2", "H3K27ac_P60_1", "H3K27ac_P0_2")



# > Differential Acetyl Analysis -------------------------------------------
# run DESeq
k27ac_dds<-DESeqDataSetFromMatrix(countData = counts, 
                                  colData = sample_sheet,
                                  design = ~ group)
k27ac_dds<-DESeq(k27ac_dds)

# format resutls
P60vP7_K27ac = results(k27ac_dds, contrast = c("group", "P60", "P7"), tidy = TRUE) %>% 
  dplyr::rename(GeneID = row) %>% 
  left_join(all_bg_sites) %>% 
  dplyr::select(-Strand) %>% 
  dplyr::mutate(k27ac_sig = case_when(log2FoldChange >= 0 & padj < 0.05 ~ "UP",
                                      log2FoldChange <= 0 & padj < 0.05 ~ "DOWN",
                                      is.na(padj) ~ "filtered",
                                      TRUE ~ "N.S."))

# save results
P60vP7_K27ac %>% write_tsv("../../results/invivo/DiffExp_H3K27ac/K27ac_data.tsv")

# > Plotting --------------------------------------------------------------
nums = table(P60vP7_K27ac$k27ac_sig)

P60vP7_K27ac %>% 
  dplyr::filter(k27ac_sig != "filtered") %>% 
  ggplot(aes(x = log10(baseMean), y = log2FoldChange, color = k27ac_sig))+
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values= c("blue", "black", "red")) +
  geom_text(x = 3, y = 5, label = paste0("Gain H3K27ac \n n = ",nums["UP"] ), color = "red")+
  geom_text(x = 3, y = -5, label = paste0("Lose H3K27ac \n n = ",nums["DOWN"] ), color = "blue")+
  theme_classic() + 
  labs( x = "log(Avg Signal)", y = "Log2FC")+
  theme(legend.position = "none")

ggsave("../../figures/H3k27ac_ma.png")  

P60vP7_K27ac %>% 
  dplyr::filter(k27ac_sig != "filtered") %>% 
  ggplot(aes(x = log10(baseMean), y = log2FoldChange, color = k27ac_sig))+
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values= c("blue", "black", "red")) +
  geom_text(x = 3, y = 5, label = paste0("Gain H3K27ac \n n = ",nums["UP"] ), color = "red")+
  geom_text(x = 3, y = -5, label = paste0("Lose H3K27ac \n n = ",nums["DOWN"] ), color = "blue")+
  theme_classic() + 
  labs( x = "log(Avg Signal)", y = "Log2FC")+
  theme(legend.position = "none") +
  xlim(c(.5, 4.3))

ggsave("../../figures/H3k27ac_ma_crop.png") 

# Gene --------------------------------------------------------------------

# > Getting gene reads ----------------------------------------------------
sample_sheet = SAMPLESHEET %>% 
  dplyr::filter(Assay_Type %in% "RNA-Seq") %>% 
  dplyr::filter(tissue %in% "whole cerebellum")

# read in gene counts for each sample 
gene.reads = list.files("../../sequencing_data/preprocessed_data/RNA/", pattern = "gene.reads", recursive = T, full.names = T) %>% 
  map_dfc(read_tsv, col_names=F)  %>% 
  # select one gene symbol column and the coulmns with the counts
  dplyr::select(c(1, seq_len(ncol(.))[ seq_len(ncol(.)) %% 2 ==0]))

colnames(gene.reads) = c("SYMBOL", list.files("../../sequencing_data/preprocessed_data/RNA/"))

# filter out the stats 
gene.reads = gene.reads[!c(grepl("__no_feature", gene.reads$SYMBOL)| 
                             grepl("__ambiguous", gene.reads$SYMBOL)| 
                             grepl("__too_low_aQual", gene.reads$SYMBOL)|  
                             grepl("__not_aligned", gene.reads$SYMBOL)| 
                             grepl("__alignment_not_unique", gene.reads$SYMBOL)), ]

head(gene.reads)
tail(gene.reads)

# > Run DESeq -------------------------------------------------------------
group = str_extract( colnames(gene.reads)[-1], "[^ _]+" )
sample_sheet$group = group
counts = gene.reads[,-1]
dds<-DESeqDataSetFromMatrix(countData = counts,
                            colData = sample_sheet,
                            design = ~ group ) 

dds<-DESeq(dds)
rownames(dds) <- gene.reads$SYMBOL

P60vP7_GE = as.data.frame(results(dds, contrast = c("group", "p60", "p7"))) %>% 
  dplyr::mutate(SYMBOL = gene.reads$SYMBOL) %>% 
  relocate(SYMBOL) %>% 
  dplyr::mutate(gene_sig = case_when(log2FoldChange >= 0 & padj < 0.05 ~ "UP",
                                      log2FoldChange <= 0 & padj < 0.05 ~ "DOWN",
                                      is.na(padj) ~ "filtered",
                                      TRUE ~ "N.S."))

# save results

as.data.frame(counts(dds, normalized = T)) %>%
  dplyr::select(starts_with(c("p7", "p60"))) %>% 
  dplyr::mutate(SYMBOL = P60vP7_GE$SYMBOL) %>% 
  dplyr::mutate(p7_mean = rowMeans(dplyr::select(., starts_with("p7"))),
         p60_mean = rowMeans(dplyr::select(.,starts_with("p60")))) %>% 
  dplyr::select(SYMBOL, p7_mean, p60_mean) %>% 
  full_join(P60vP7_GE) %>% write_tsv("../../results/invivo/DiffExp_RNA/GeneExp_data.tsv")


as.data.frame(counts(dds, normalized = T)) %>% dplyr::mutate(SYMBOL = rownames(.)) %>%  relocate(SYMBOL) %>% write_tsv("../../results/invivo/DiffExp_RNA/normalized_counts.tsv")
as.data.frame(counts(dds, normalized = F)) %>% dplyr::mutate(SYMBOL = rownames(.)) %>%  relocate(SYMBOL) %>% write_tsv("../../results/invivo/DiffExp_RNA/raw_counts.tsv")

# > Plotting --------------------------------------------------------------
nums = table(P60vP7_GE$gene_sig)

P60vP7_GE %>% 
  dplyr::filter(gene_sig != "filtered") %>% 
  ggplot(aes(x = log10(baseMean), y = log2FoldChange, color = gene_sig))+
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values= c("blue", "black", "red")) +
  geom_text(x = 3, y = 5, label = paste0("Up-regulated \n n = ",nums["UP"] ), color = "red")+
  geom_text(x = 3, y = -5, label = paste0("Down-regulated \n n = ",nums["DOWN"] ), color = "blue")+
  theme_classic() + 
  labs( x = "log(Avg Signal)", y = "Log2FC")+
  theme(legend.position = "none")

ggsave("../../figures/GE_ma.png")  

P60vP7_GE %>% 
  dplyr::filter(!is.na(padj)) %>% 
  dplyr::mutate(sig = case_when(abs(log2FoldChange) >= 1 & padj < 0.05 ~ "FDR < 0.05 & |Log2FC| >= 1",
                                abs(log2FoldChange) >=1 & padj > 0.05 ~ "|Log2FC| >=1",
                                abs(log2FoldChange) < 1 & padj < 0.05 ~ "FDR < 0.05",
                                abs(log2FoldChange) < 1 & padj > 0.05 ~ "N.S.")) %>% 
  dplyr::mutate(label = case_when( padj < 5e-50 ~ SYMBOL)) %>% 
  ggplot(aes(y = -log10(padj), x = log2FoldChange, color = sig, label = label))+
  geom_point( alpha = 0.5) +
  scale_color_manual(values= c("black", "dodgerblue2", "deeppink4", "grey")) +
  geom_text( check_overlap = T, color = "black", vjust = 0, hjust = 1, fontface = "bold") +
  geom_vline( xintercept = c(-1, 1), color = "grey", linetype="dashed")+
  geom_hline( yintercept = -log10(0.05), color = "grey", linetype="dashed")+
  theme_classic() + 
  labs( x = "log2FC", y = "-log(adj. p-value)")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA))


ggsave("../../figures/GE_volcano.png")  




