# Loading packages --------------------------------------------------------
library(tidyverse)
library("ggpmisc")
library(GenomicRanges)
library(rtracklayer)
library(Gviz) 
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(biomaRt)
library(GenomicInteractions)
options(ucscChromosomeNames=FALSE) 


# Read in data ------------------------------------------------------------

# > Chromatin Links   ------------------------------------------------------
# Source: Yamada et al. 2019 https://www.nature.com/articles/s41586-019-1190-7
plac_file <- read_delim("../../sequencing_data/Yamada/combined_MAPS_peaks.txt", 
                        delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)



# > Peaksets   -------------------------------------------------------------
DNaseP7 <- read_delim("../../results/mergedPeaks/DNaseP7.bed", 
                      delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

DNaseP60 <- read_delim("../../results/mergedPeaks/DNaseP60.bed", 
                       delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

H3K27acP7 <- read_delim("../../results/mergedPeaks/H3K27acP7.bed", 
                        delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

H3K27acP60 <- read_delim("../../results/mergedPeaks/H3K27acP60.bed", 
                         delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

ZicP7 <- read_delim("../../results/mergedPeaks/ZicP7.bed", 
                    delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

ZicP60 <- read_delim("../../results/mergedPeaks/ZicP60.bed", 
                     delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# Define Functions --------------------------------------------------------
plot_track <- function(gene, expand_l= 0, expand_r = 0, ylim = c(0,3), save = F){
  # getting gene info
  id = mapIds(org.Mm.eg.db, keys =gene , column ="ENTREZID", keytype="SYMBOL")
  info = select(TxDb.Mmusculus.UCSC.mm10.knownGene, keys = id  , columns=c("TXNAME", "TXSTART", "TXEND", "TXCHROM"), keytype="GENEID")
  start = info$TXSTART[1]
  end = info$TXEND[1]
  chr = info$TXCHROM[1]
  
  # setting up plot inputs
  start = start - expand_l
  end = end + expand_r  
  
  
  tracks <- list(interaction_track,
                 p7_zic_track,
                 p60_zic_track,
                 p7_dnase_track,
                 #p60_dnase_track,
                 p7_h3k27ac_track,
                 p60_h3k27ac_track,
                 p7_dt,
                 p7_2dt,
                 p60_dt,
                 p60_2dt,
                 customFromTxDb,
                 genomeAxis)
  
  #plot
  
 
  
  if(save){
    png(paste0("../../figures/",gene,"_loops.png"), units = "in", res = 300, height = 2, width = 6)
    
    plotTracks(tracks,  
               from=start,
               to=end, 
               chromosome=chr,
               type="hist",
               window = 1000,
               ylim = ylim,
               sizes = rep(1, length(tracks)), 
               transcriptAnnotation = "symbol",
               rotation.title = 1)
    
    dev.off()
  }else{
    plotTracks(tracks,  
               from=start,
               to=end, 
               chromosome=chr,
               type="hist",
               window = 1000,
               ylim = ylim,
               sizes = rep(1, length(tracks)), 
               transcriptAnnotation = "symbol",
               rotation.title = 1)
  }
  
  
  
}




# Wrangle data ------------------------------------------------------------

# prepare data
p7_bw <- import.bw("../../sequencing_data/preprocessed_data/zic_chip/SRR1557091/SRR1557091_norm.bw", as="GRanges") 
p7_2bw <- import.bw("../../sequencing_data/preprocessed_data/zic_chip/SRR1557092/SRR1557092_norm.bw", as="GRanges") 

p60_bw <- import.bw("../../sequencing_data/preprocessed_data/zic_chip/SRR1557093/SRR1557093_norm.bw", as="GRanges") 
p60_2bw <- import.bw("../../sequencing_data/preprocessed_data/zic_chip/SRR1557094/SRR1557094_norm.bw", as="GRanges") 

genomeAxis <- GenomeAxisTrack(name="MyAxis") 

customFromTxDb <- GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene) 

p7_dt <- DataTrack(p7_bw, genome = "mm10",name = "Zic P7") 
p7_2dt <- DataTrack(p7_2bw, genome = "mm10",name = "Zic P7") 

p60_dt <- DataTrack(p60_bw, genome = "mm10", name = "Zic P60") 
p60_2dt <- DataTrack(p60_2bw, genome = "mm10", name = "Zic P60") 

p7_zic_track = AnnotationTrack(GRanges(seqnames = ZicP7$X1, ranges = IRanges(ZicP7$X2, ZicP7$X3)), 
                               genome = "mm10", name = "P7 Zic", fill = "blue")
p60_zic_track = AnnotationTrack(GRanges(seqnames = ZicP60$X1, ranges = IRanges(ZicP60$X2, ZicP60$X3)),
                                genome = "mm10", name = "P60 Zic", fill = "red")

p7_dnase_track = AnnotationTrack(GRanges(seqnames = DNaseP7$X1, ranges = IRanges(DNaseP7$X2, DNaseP7$X3)), 
                                 genome = "mm10", name = "P7n Dnase", fill = "blue")
p60_dnase_track = AnnotationTrack(GRanges(seqnames = DNaseP60$X1, ranges = IRanges(DNaseP60$X2, DNaseP60$X3)), 
                                  genome = "mm10", name = "P60 Dnase", fill ="red")

p7_h3k27ac_track = AnnotationTrack(GRanges(seqnames = H3K27acP7$X1, ranges = IRanges(H3K27acP7$X2, H3K27acP7$X3)), 
                                   genome = "mm10", name = "P7 H3K27ac", fill = "blue")
p60_h3k27ac_track = AnnotationTrack(GRanges(seqnames = H3K27acP60$X1, ranges = IRanges(H3K27acP60$X2, H3K27acP60$X3)), 
                                    genome = "mm10", name = "P60 H3k27ac", fill ="red")

interaction_track <- InteractionTrack(GenomicInteractions(GRanges(seqnames = plac_file$X1, ranges = IRanges(plac_file$X2, plac_file$X3)), 
                                                          GRanges(seqnames = plac_file$X4, ranges = IRanges(plac_file$X5, plac_file$X6))), 
                                      name = "H3K4me3 PLAC-seq")


# Plot --------------------------------------------------------------------

plot_track("Grin2c", expand_l = 100000, expand_r = 1000)

plot_track("Wnt7b",  expand_l = 100000, expand_r = 10000)






