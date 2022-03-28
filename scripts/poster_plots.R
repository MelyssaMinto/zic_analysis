library(tidyverse)
library(cowplot)


# Read in data ------------------------------------------------------------



# figure1 -----------------------------------------------------------------
gene_expr = read_tsv("../results/DiffExp_RNA/normalized_counts.tsv")
zic1_kd =    read_excel("../sequencing_data/Frank_2015/41593_2015_BFnn3995_MOESM88_ESM.xls", sheet = 1)
zic2_kd =    read_excel("../sequencing_data/Frank_2015/41593_2015_BFnn3995_MOESM88_ESM.xls", sheet = 2)
zic_p60vp7 = read_tsv("../results/DiffExp_ZicChIP/ZicChIPDA_data.tsv")
zic_3v7 =  read_excel("../sequencing_data/Frank_2015/41593_2015_BFnn3995_MOESM82_ESM.xls", sheet = "DIV3_v_DIV7_genes.txt")

#> panel a -----------------------------------------------------------------


a = rasterGrob(readPNG("../figures/cerebellum_dev.png"), interpolate = T)


#> panel b -----------------------------------------------------------------



b = gene_expr %>% 
  dplyr::filter( grepl("Zic", SYMBOL)) %>% 
  pivot_longer(cols = starts_with("p")) %>% 
  dplyr::mutate(time = str_extract(name, "[^\\_]+")) %>% 
  ggplot(aes(x = time, y = value, fill = time)) +
  scale_fill_manual(values = c("blue", "orange", "red"),limits = c("p7","p14", "p60") )+
  scale_x_discrete( labels = c("P7", "P14", "P60"), limits = c("p7","p14", "p60")) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 50000))+
  stat_summary(fun.data="mean_se",geom="errorbar", color = "black", size = .5, width = .5) +
  stat_summary(fun="mean",geom="bar") +
  facet_wrap(~SYMBOL, ncol = 5) +
  labs(x = "", y = "normalized counts")+
  theme_bw() +
  theme( panel.grid = element_blank(),
         legend.position = "none") 



#> panel c -----------------------------------------------------------------


D3vD7 <-zic_3v7 %>% 
  dplyr::select(gene:q_value) %>% 
  dplyr::mutate(log2.fold_change. = as.numeric(log2.fold_change. ),
         test_stat = as.numeric(test_stat ) ) %>% 
  #dplyr::filter(q_value < 0.05) %>% 
  dplyr::mutate(LFC_WT = log2.fold_change.) %>% 
  dplyr::select(gene, LFC_WT)  

Zic1KD_D7 <- zic1_kd %>% 
  #plyr::filter(q_value < 0.05) %>% 
  dplyr::mutate(LFC_Zic1KD = `log2(Fold-Change)`) %>% 
  dplyr::select(gene, LFC_Zic1KD) 

Zic2KD_D7 <- zic2_kd %>% 
  #dplyr::filter(q_value < 0.05) %>% 
  dplyr::mutate(LFC_Zic2KD = `log2(Fold-Change)`) %>% 
  dplyr::select(gene, LFC_Zic2KD)


#> combining data 


c = D3vD7 %>% 
  full_join(Zic1KD_D7) %>% 
  drop_na() %>% 
  ggplot(aes(x = LFC_WT, y = LFC_Zic1KD, label = gene)) +
  geom_point() +
  # geom_smooth(method = "lm", se=FALSE)+
  # stat_regline_equation(label.y = 4.5, aes(label = ..eq.label..), color = "black") +
  # stat_regline_equation(label.y = 4, aes(label = ..rr.label..), color = "black") +
  # geom_text_repel()+
  labs(x ="LFC DIV7/DIV3", y = "LFC Zic1/ +7 DIV Control") +
  theme_bw() +
  ylim(-6, 6) +
  xlim(-6, 6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank())

D3vD7 %>% 
  full_join(Zic2KD_D7) %>% 
  drop_na() %>% 
  ggplot(aes(x = LFC_WT, y = LFC_Zic2KD, label = gene)) +
  geom_point() +
  # geom_smooth(method = "lm", se=FALSE)+
  # stat_regline_equation(label.y = 4.5, aes(label = ..eq.label..), color = "black") +
  # stat_regline_equation(label.y = 4, aes(label = ..rr.label..), color = "black") +
  # geom_text_repel()+
  labs(x ="LFC DIV7/DIV3", y = "LFC Zic2/ +7 DIV Control") +
  theme_bw() +
  ylim(-6, 6) +
  xlim(-6, 6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.grid = element_blank())





#> panel d -----------------------------------------------------------------

nums = table(zic_p60vp7$zic_sig)
d = zic_p60vp7 %>% 
  dplyr::filter(zic_sig != "filtered") %>% 
  ggplot(aes(x = log10(baseMean), y = log2FoldChange, color = zic_sig))+
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values= c("blue", "black", "red")) +
  geom_text(x = 3.5, y = 5, label = paste0("Increased Zic Binding\n n = ",nums["UP"] ), color = "red")+
  geom_text(x = 3.5, y = -6, label = paste0("Decreased Zic Binding\n n = ",nums["DOWN"] ), color = "blue")+
  theme_classic() + 
  xlim(c(.7, 4.7))+
  labs( x = "log10(Avg Zic ChIP  Signal)", y = "Log2FC")+
  theme(legend.position = "none") 



#> Figure ------------------------------------------------------------------
png("../figures/poster_fig_1.png", height = 8, width = 10, res = 300, units = "in")
plot_grid(a,b,c,d, nrow = 2, labels = LETTERS[1:4])
dev.off()


# Figure2 -----------------------------------------------------------------

P60vP7_zic = read_tsv("../results/DiffExp_ZicChIP/ZicChIPDA_data.tsv")

workflow = readPNG('../figures/homer_woorkflow.png')

homer_results_P7 = read_tsv("../results/homer_results/P60vP7_DOWN/homer_results.tsv")
homer_results_P60 = read_tsv("../results/homer_results/P60vP7_UP/homer_results.tsv")
homer_results_static = read_tsv("../results/homer_results/P60vP7_NS/homer_results.tsv")


# > panel a ---------------------------------------------------------------

a = rasterGrob(workflow, interpolate = T)

# >panel b ----------------------------------------------------------------

b =  bind_rows(
  list(
    #   # static motifs that are transcriptionally enriched
    # static = homer_results_static %>%
    # dplyr::filter(`q-value (Benjamini)` < 0.05) %>%
    # dplyr::filter(padj < 0.05 ) %>% 
    # dplyr::filter(abs(log2FoldChange) <= 1 ) %>% 
    # dplyr::filter(name != "Unknown-ESC-element") %>%
    # dplyr::mutate(label = paste0(name, "(", SYMBOL, ")")),
    # p60 motifs that are transcriptionally enriched
    p60 = homer_results_P60 %>%
      dplyr::filter(`q-value (Benjamini)` < 0.05) %>%
      dplyr::filter(padj < 0.05 ) %>% 
      dplyr::filter(log2FoldChange >= 1 ) %>% 
      dplyr::filter(name != "Unknown-ESC-element") %>%
      dplyr::mutate(label = ifelse(str_to_lower(name) != str_to_lower(SYMBOL), paste0(name, "(", SYMBOL, ")"), name)),
    # p7 motifs that are transcriptionally enriched
    p7 = homer_results_P7 %>%
      dplyr::filter(`q-value (Benjamini)` < 0.05) %>%
      dplyr::filter(padj < 0.05 ) %>% 
      dplyr::filter(log2FoldChange <= -1 ) %>% 
      dplyr::filter(name != "Unknown-ESC-element") %>%
      dplyr::mutate(label = ifelse(str_to_lower(name) != str_to_lower(SYMBOL), paste0(name, "(", SYMBOL, ")"), name))), 
  .id = "id")  %>% 
  dplyr::mutate(id = factor(id, levels = c("p7", "static", "p60"), labels = c('P7 Zic Peaks','Static Zic Peaks','P60 Zic peaks' ))) %>% 
  dplyr::arrange(desc(id), `q-value (Benjamini)`) %>% 
  dplyr::mutate(label_f = factor(label, unique(label))) %>% 
  dplyr::filter(baseMean >= 100) %>% 
  ggplot(aes(y = log2FoldChange, x = label_f, fill = `q-value (Benjamini)`)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0,color = "grey",linetype = "dashed") +
  #coord_flip() +
  labs(x = "TF (mapped gene)", fill = "Motif Enrichment\n q-value") +
  scale_y_continuous(limits = c(-11,5), expand = c(0, 0))+
  facet_grid(.~id, scales = "free_x", space = "free_x")+
  theme_classic() +
  theme( axis.line.y = element_blank(),
         axis.text.x = element_text(angle = 45, hjust = 1, size =12),
         legend.position = "none") +
  ylab(bquote("P60/P7" ~ log[2] ~ "Fold Enrichment"))



# > panel c/d -------------------------------------------------------------
plot_fam <- function(data, xlab, regulation, baseMeanCutoff = 100){
  
  
  data %>% 
    dplyr::filter(baseMean >= baseMeanCutoff) %>% 
    dplyr::filter(`q-value (Benjamini)`<0.05) %>% 
    dplyr::mutate(reg = case_when(log2FoldChange >= 1 ~ "UP", 
                                  log2FoldChange <= -1 ~ "DOWN",
                                  TRUE ~ "N.S.") ) %>% 
    dplyr::filter(reg %in% regulation) %>% 
    dplyr::select(fam) %>% 
    # cleaning up family names
    dplyr::mutate(fam = strsplit(as.character(fam), ",")) %>% 
    unnest(fam) %>% 
    dplyr::mutate(fam = strsplit(as.character(fam), ":")) %>% 
    unnest(fam) %>% 
    dplyr::filter(!is.na(fam)) %>% 
    dplyr::filter(!(fam %in% c("?"))) %>% 
    group_by(fam) %>% 
    dplyr::count() %>% 
    ungroup() %>% 
    dplyr::arrange(n) %>% 
    dplyr::mutate(fam=factor(fam, levels=fam)) %>% 
    ggplot(aes(x = fam, y = n))+
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_classic() +
    scale_y_continuous( expand = c(0, 0)) +
    labs( y = xlab, x = "" )
} 

c = plot_fam(homer_results_P7, xlab = "# of TF Families enriched at P7", regulation = "DOWN", baseMeanCutoff = 100)
d = plot_fam(homer_results_P60,  xlab = "# of TF Families enriched at P60", regulation = "UP", baseMeanCutoff = 100)


# > figure ----------------------------------------------------------------
png("../figures/poster_fig_2.png", height =6, width = 10, res = 300, units = "in")
plot_grid(plot_grid(a,b, ncol = 1, labels = LETTERS[1:2], rel_heights = c(.3, .7)),
          plot_grid(c,d, ncol = 1 , labels = LETTERS[3:4]),
          ncol = 2,
          rel_widths = c(.7, .3)
)
dev.off()


# Figure 3 ----------------------------------------------------------------

atoh1_zic_overlap_p7 = read_tsv("../results/mergedPeaks/zic_atoh1_p7.bed", col_names = F) 
atoh1_zic_overlap_p60 = read_tsv("../results/mergedPeaks/zic_atoh1_p60.bed", col_names = F) 
atoh1_zic_overlap_static = read_tsv("../results/mergedPeaks/zic_atoh1_static.bed", col_names = F) 

atoh_peaks =  read_tsv("../results/mergedPeaks/atoh1.bed", col_names = F)

P60vP7_zic_overlap = read_tsv("../results/peak_gene/zic_atoh_overlap.txt")

zic_p7_overlap = atoh1_zic_overlap_p7 %>% 
  dplyr::select(X1, X2, X3) %>% 
  distinct()
atoh1_p7_overlap = atoh1_zic_overlap_p7 %>% 
  dplyr::select(X4, X5, X6) %>% 
  distinct()

zic_p60_overlap = atoh1_zic_overlap_p60 %>% 
  dplyr::select(X1, X2, X3) %>% 
  distinct()
atoh1_p60_overlap = atoh1_zic_overlap_p60 %>% 
  dplyr::select(X4, X5, X6) %>% 
  distinct()


zic_static_overlap = atoh1_zic_overlap_static %>% 
  dplyr::select(X1, X2, X3) %>% 
  distinct()
atoh1_static_overlap = atoh1_zic_overlap_static %>% 
  dplyr::select(X4, X5, X6) %>% 
  distinct()

# > panel a ---------------------------------------------------------------
a = data.frame(group = c("P7", "P60", "Static"),
           count = c(nrow(atoh1_p7_overlap),nrow(atoh1_p60_overlap), nrow(atoh1_static_overlap) )/ nrow(atoh_peaks)) %>% 
  add_row(group = "No overlap", count = 1 - sum(c(nrow(atoh1_p7_overlap),nrow(atoh1_p60_overlap), nrow(atoh1_static_overlap) )/ nrow(atoh_peaks))) %>% 
  dplyr::mutate(col = "Atoh1 Peak") %>% 
  ggplot(aes(x = col, y = count, fill = group)) +
  geom_col() +
  scale_fill_manual(limits = c("P7", "P60", "Static", "No Overlap"), labels = c("P7", "P60", "Static", "No \nOverlap"),values = c("blue", "red", "black", "grey")) +
  theme_classic() +
  labs(x = "", y = "% Overlap0", fill = "Zic Peak") 

# > panel b ---------------------------------------------------------------
# What proportion of the Zic peaks overlap with the atoh1 pwaks
# create data frame
mat  = data.frame(overlap = c(nrow(zic_p7_overlap),nrow(zic_p60_overlap), nrow(zic_static_overlap) ),
                  npeaks = c(sum(P60vP7_zic$zic_sig == "DOWN"), sum(P60vP7_zic$zic_sig == "UP"), sum(P60vP7_zic$zic_sig == "N.S."))) %>%
  dplyr::mutate( non_overlap = npeaks - overlap) %>% 
  dplyr::select(-npeaks)

row.names(mat) =  c("P7", "P60", "Static")

# compute stats
p7vp60  = prop.test(as.matrix(mat)[1:2,])
p7vstatic  = prop.test(as.matrix(mat)[c(1,3),])
p60vstatic  = prop.test(as.matrix(mat)[c(2,3),])

#plot
b = mat %>% 
  dplyr::mutate(group = rownames(.),
                group = fct_relevel(group, c("P7", "Static", "P60"))) %>% 
  pivot_longer(cols = c("overlap", "non_overlap")) %>% 
  ggplot(aes(x = group, y = value, fill = name)) +
  geom_col(position = "fill") +
  labs(x = "Zic Peak", y = "% Overlap", fill = "Overlap\n w/Atoh1") +
  scale_fill_manual(limits = c("non_overlap", "overlap"), labels = c("No", "Yes"), values = c("black", "grey")) +
  geom_signif(tip_length = .000001, annotations = c("***", "***"), xmin = c(0.8, 0.8), xmax = c(2.2,3.2), y_position = c(1.01, 1.1) , color = "black") +
  theme_classic() +
  scale_y_continuous(expand = expansion(add = c(0, .2)), limits = c(0,1), oob = rescale_none)
b
        
# > panel c ---------------------------------------------------------------
cont = P60vP7_zic_overlap %>% 
  dplyr::filter(Overlap %in% c(TRUE)) %>% 
  dplyr::filter(zic_sig %in% c("UP", "DOWN", "N.S.")) %>% 
  dplyr::select(zic_sig, group)

tab = table(cont$zic_sig, cont$group)


sig_all = chisq_test(tab)
sig_p7vstatic = chisq_test(tab[c("DOWN", "N.S."),])
sig_p7vp60 = chisq_test(tab[c("DOWN", "UP"),])
sig_staticvp60 = chisq_test(tab[c("N.S.", "UP"),])

c = P60vP7_zic_overlap %>% 
  dplyr::filter(Overlap %in% c(TRUE)) %>% 
  dplyr::filter(zic_sig %in% c("UP", "DOWN", "N.S.")) %>% 
  dplyr::count(zic_sig, group) %>% 
  ggplot(aes( x = zic_sig, y = n,  fill = group)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "% Peaks", x = "Overlap with Atoh1", fill = "Genomic Region")+
  scale_x_discrete(limits = c("DOWN", "N.S.", "UP"), labels = c("P7", "Static", "P60"))+
  theme_classic() +
  ggsignif::geom_signif(tip_length = .0000001, annotations = sig_all$p.signif, xmin = c(0.8), xmax = c(3.2), y_position = c(1.1) , color = "black")
c
# > panel d ---------------------------------------------------------------

d2 = rasterGrob(readPNG("../figures/fig1_track2.png"))
d1 = rasterGrob(readPNG("../figures/fig1_track1.png"))


# > figure ----------------------------------------------------------------
png("../figures/poster_fig_3.png", height = 8, width = 10, res = 300, units = "in")
plot_grid(plot_grid(a,b,c, ncol = 3, labels = LETTERS[1:3], rel_widths = c(.25, .3, .45)),
          plot_grid(d1, d2, ncol = 1, labels = c("D", ""), vjust = -1),
          ncol = 1,
          #scale = c(1, 1.1), 
          rel_heights = c(1, 2))

dev.off()
# Figure 4 ----------------------------------------------------------------
# > data ------------------------------------------------------------------


# Mapped peaks to genes via chromatin loops/links
mapped_data = read_tsv("../results/FinalTables/mapped_data.txt")

# > panel a ---------------------------------------------------------------
a = rasterGrob(readPNG("../figures/fig2_d.png"))

# > panel b ---------------------------------------------------------------
# how many anchors are mapped to genes
b = mapped_data %>% 
  dplyr::select(loop_id, gene_name, gene_sig) %>% 
  dplyr::filter(gene_sig != "filtered") %>%
  dplyr::mutate(in_anchor = ifelse(is.na(loop_id), FALSE, TRUE)) %>% 
  dplyr::select(-loop_id) %>% 
  distinct() %>% 
  group_by(in_anchor) %>% 
  dplyr::count(gene_sig) %>% 
  ggplot(aes(x = gene_sig, y = n, label = n))+
  geom_col() +
  geom_text(vjust = 0, fontface = "bold", color = "black", size = 3)+
  facet_wrap(~in_anchor, labeller = labeller(in_anchor = c("FALSE" = "Not in Anchor", "TRUE" = "With in Anchor"))) +
  labs(y = "# genes \nmapped to anchors", x = "Gene Regulation") +
  theme_bw()+
  theme(plot.title = element_text(face = "bold"))

# > panel c ---------------------------------------------------------------
# How many zic peaks are not in anchors
c = mapped_data %>% 
  dplyr::select(loop_id, zic_peak, zic_sig) %>% 
  dplyr::filter(zic_sig != "filtered") %>%
  dplyr::filter(!grepl("chr[X|Y]", zic_peak)) %>% 
  dplyr::mutate(in_anchor = ifelse(is.na(loop_id), FALSE, TRUE)) %>% 
  dplyr::select(-loop_id) %>% 
  distinct() %>% 
  group_by(in_anchor) %>% 
  dplyr::count(zic_sig) %>% 
  ggplot(aes(x = zic_sig, y = n, label = n))+
  geom_col() +
  geom_text(vjust = 0, fontface = "bold", color = "black", size = 3)+
  facet_wrap(~in_anchor, labeller = labeller(in_anchor = c("FALSE" = "Not in Anchor", "TRUE" = "With in Anchor"))) +
  labs(y = "# Zic peaks \nmapped to anchors", x = "Zic Regulation") +
  theme_bw()+
  theme(plot.title = element_text(face = "bold"))
# > panel d ---------------------------------------------------------------
# How many DNase peaks are not in anchors
d = mapped_data %>% 
  dplyr::select(loop_id, dnase_peak, dnase_sig) %>%
  dplyr::filter(dnase_sig != "filtered") %>% 
  dplyr::filter(!grepl("chr[X|Y]", dnase_peak)) %>% 
  dplyr::mutate(in_anchor = ifelse(is.na(loop_id), FALSE, TRUE)) %>% 
  dplyr::select(-loop_id) %>% 
  distinct() %>% 
  group_by(in_anchor) %>% 
  dplyr::count(dnase_sig) %>% 
  ggplot(aes(x = dnase_sig, y = n, label = n))+
  geom_col() +
  geom_text(vjust = 0, fontface = "bold", color = "black", size = 3)+
  facet_wrap(~in_anchor, labeller = labeller(in_anchor = c("FALSE" = "Not in Anchor", "TRUE" = "With in Anchor"))) +
  labs(y = "# DHS peaks \nmapped to anchors", x = "DNase Regulation") +
  theme_bw()+
  theme(plot.title = element_text(face = "bold"))
# > panel e ---------------------------------------------------------------

# How many H3K27ac peaks are not in anchors
e = mapped_data %>% 
  dplyr::select(loop_id, k27ac_peak, k27ac_sig) %>%
  dplyr::filter(k27ac_sig != "filtered") %>% 
  dplyr::filter(!grepl("chr[X|Y]", k27ac_peak)) %>% 
  dplyr::mutate(in_anchor = ifelse(is.na(loop_id), FALSE, TRUE)) %>% 
  dplyr::select(-loop_id) %>% 
  distinct() %>% 
  group_by(in_anchor) %>% 
  dplyr::count(k27ac_sig) %>% 
  ggplot(aes(x = k27ac_sig, y = n, label = n))+
  geom_col() +
  geom_text(vjust = 0, fontface = "bold", color = "black", size = 3)+
  facet_wrap(~in_anchor, labeller = labeller(in_anchor = c("FALSE" = "Not in Anchor", "TRUE" = "With in Anchor"))) +
  labs(y = "# H3K27ac peaks \nmapped to anchors", x = "H3K27ac Regulation") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))




# > fig -------------------------------------------------------------------
png("../figures/poster_fig_4.png", height = 6, width = 8, res = 300, units = "in")

plot_grid(a, 
          plot_grid(b,c,d,e,nrow = 2, labels = LETTERS[2:5]),
        
          labels = c("A", ""),
          ncol = 1, 
          vjust = -.5)


dev.off()




# Figure 5 ----------------------------------------------------------------
# > panel a ---------------------------------------------------------------
a = mapped_data %>% 
  drop_na(c("loop_id", "zic_peak", "gene_sig", "zic_sig")) %>% 
  dplyr::select(zic_peak, gene_name, zic_sig, gene_sig) %>% 
  distinct() %>% 
  dplyr::filter( zic_sig != "filtered") %>% 
  dplyr::filter( gene_sig != "filtered") %>% 
  dplyr::count(zic_sig, gene_sig) %>% 
  ggplot(aes(x = zic_sig, y =n,  fill = gene_sig))+
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label=n),position = position_fill(vjust = 0.5), color = "white", fontface = "bold", size = 6) +
  scale_x_discrete(labels=c("DOWN" = "P7 Peak", "N.S." = "Static", "UP" = "P60 Peak"), expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual( labels = c("Constituitive", "Down-regulated", "Up-regulated"), values = c("black", "blue", "red"))+
  theme_classic() +
  labs(x = "Differential Signal of Zic", fill = "Diff Signal", y = "% Genes") +
  theme( legend.position = c(.8,.09),
         legend.text = element_text(size = 5),
         legend.title = element_text(size = 8, face = "bold"),
         legend.key.size = unit(.3, "lines"),
         legend.box.background = element_rect(color = "black"))
         
# > panel b ---------------------------------------------------------------
df = mapped_data %>%
  dplyr::select(gene_name, zic_sig, p60_mean, p7_mean)  %>%
  dplyr::filter(zic_sig != "filtered") %>%
  group_by(zic_sig) %>%
  distinct() %>%
  dplyr::filter(!is.na(zic_sig)) %>%
  pivot_longer(cols = ends_with("mean"), values_to = "means", names_to = "timepoint") %>%
  dplyr::mutate(timepoint = factor(gsub("_mean", "", timepoint)))

stat.test = df %>%
  t_test(means ~ timepoint) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(fun = "mean_se")


b = df %>% 
  group_by(zic_sig, timepoint) %>% 
  summarise( mean= mean(means, na.rm = T), 
             sd = sd(means, na.rm = T)/sqrt(n()),
             .groups = "drop") %>% 
  ungroup() %>% 
  # plotting
  dplyr::mutate(timepoint = fct_relevel(timepoint, levels = c("p7", "p60") )) %>% 
  ggplot(aes(x = zic_sig, y = mean, ymin = mean - sd, ymax = mean + sd, fill = timepoint, color = timepoint)) +
  geom_errorbar(width = 0.2, position = position_dodge(.8)) +
  scale_fill_manual(values = c( "blue", "red"),labels = c("p60" = "P60", "p7" = "P7"))+
  scale_color_manual(values = c( "blue", "red"),labels = c("p60" = "P60", "p7" = "P7"))+
  scale_x_discrete(labels=c("DOWN" = "P7 Peak", "N.S." = "Static", "UP" = "P60 Peak"), expand = c(0,0))+
  theme_classic()+
  xlab("Zic Peak Enrichement") +
  ylab("Mean Expresson of Mapped Gene") +
  ggsignif::geom_signif(
    tip_length = 0,
    annotations = stat.test$p.adj.signif,
    xmin = c(0.8, 1.8, 2.8),
    xmax = c(1.2,2.2,3.2), 
    y_position = stat.test$y.position,
    color = "black") +
  geom_point(size = 1, position = position_dodge(0.8)) +
  theme(legend.position="bottom") +
  labs(fill = "Expression", color = "Expression")

# > panel c ---------------------------------------------------------------
n_peaks_to_gene = mapped_data %>% 
  drop_na("loop_id") %>% 
  dplyr::select(zic_peak, gene_name) %>% 
  distinct() %>% 
  dplyr::filter(!is.na(zic_peak)) %>% 
  group_by(gene_name) %>% 
  dplyr::count() %>% 
  ungroup() %>% 
  dplyr::filter(!is.na(gene_name)) %>% 
  left_join(mapped_data %>% dplyr::select(gene_name, gene_sig)) %>% 
  dplyr::mutate(gene_name = fct_reorder(gene_name, n)) %>% 
  dplyr::filter(gene_sig != "filtered") %>% 
  distinct() 

top_genes = n_peaks_to_gene %>% 
  dplyr::filter(gene_sig != "N.S.") %>% 
  dplyr::arrange(n) %>% 
  group_by(gene_sig) %>% 
  top_n(70, n) %>% 
  pull(gene_name)


c = mapped_data %>% 
  drop_na("loop_id") %>% 
  dplyr::select(zic_peak, gene_name) %>% 
  distinct() %>% 
  dplyr::filter(!is.na(zic_peak)) %>% 
  group_by(gene_name) %>% 
  dplyr::add_count() %>% 
  ungroup() %>% 
  left_join(mapped_data %>% dplyr::select(gene_name, gene_sig, zic_peak, zic_sig) %>%  distinct()) %>% 
  dplyr::filter(gene_name %in% top_genes) %>% 
  dplyr::mutate(gene_name = fct_reorder(gene_name, n))  %>% 
  dplyr::arrange(desc(n)) %>%
  ggplot(aes(y = gene_name, fill =zic_sig))+
  geom_bar(stat = "count", ) +
  #coord_flip()+
  # geom_text(aes(label = n), hjust = 1, color = "white", fontface = "bold") +
  facet_wrap(.~gene_sig, scales = "free", labeller = labeller(gene_sig = c("DOWN" = "Down-regulated", "UP" = "Up-regulated")), ncol = 2)+
  scale_fill_manual( name = "Zic Diff Signal", labels = c("P7 Peak", "Static", "P60 Peak"), values = c("blue", "black", "red"))+
  theme_classic() + 
  labs(y = "Gene Symbol", x = "# Zic Peaks Mapped to Gene") +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1), legend.position = c(.99, .25),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.box.background = element_rect(color="black", size=1))



c

#  >fig -------------------------------------------------------------------
png("../figures/poster_fig_5.png", height = 8, width = 10, res = 300, units = "in")
plot_grid(plot_grid(a,b, ncol = 1, labels = c("A", "B")),
          c,
          ncol = 2, 
          labels = c("", "C"),
          rel_widths = c(.3, .7))

dev.off()

# Figure 6 ----------------------------------------------------------------

# >data -------------------------------------------------------------------


earlyVlate_act = read_tsv("../results/peak_gene/rrho/early_late_activating.txt")
early_actVrep = read_tsv("../results/peak_gene/rrho/early_activating_repressive.txt")
earlyVlate_rep = read_tsv("../results/peak_gene/rrho/early_late_repressive.txt")
late_actvrep = read_tsv("../results/peak_gene/rrho/late_activating_repressive.txt")
actvrep = read_tsv("../results/peak_gene/rrho/activating_repressive.txt")
earlyvlate = read_tsv("../results/peak_gene/rrho/early_late.txt")
early_allVloop = read_tsv("../results/peak_gene/rrho/early_allvloop.txt")
late_allVloop = read_tsv("../results/peak_gene/rrho/late_allvloop.txt")

plot_distinct <- function(data, cond1, cond2, label1, label2) {
  dis = data %>% 
    dplyr::filter(gene_sig %in% c("UP", "DOWN")) %>% 
    
    dplyr::mutate(label = case_when(enriched_in %in% "Similar Enrichment" ~ "" ,
                                    TRUE ~ str_to_title(TF))) %>% 
    dplyr::select(-SYMBOL, -gene_sig, -ends_with("mean")) %>% 
    distinct() 
  
  dis_table = dis$enriched_in %>% 
    fct_relabel(~ str_to_title(gsub("_", " ", .x))) %>% 
    fct_count()
  colnames(dis_table) = c("Enrichment", "n")
  
  dis %>% 
    ggplot(aes(x = -log10(get(cond1)), y = -log10(get(cond2)), color = enriched_in, label = label )) +
    geom_point(size = 0.5) +
    geom_text_repel(show.legend = F, max.overlaps = 20) +
    theme_classic() +
    labs(x = paste0(gsub("_", " ", cond1), " \n-log10(adj p-val)"), y = paste0(gsub("_", " ", cond2), " \n-log10(adj p-val)") ) +
    scale_color_manual(name = "", limits = c(cond1, cond2, "Similar Enrichment"), labels = c(label1, label2, "Similar \nEnrichment"), values = c("dodgerblue4", "darkorange4", "grey") ) +
    annotate(geom = "table", label = list(dis_table), x = 4.5, y = 4.5, fill = "white" ,  table.theme = ttheme_gtlight)+
    scale_y_continuous(limits = c(0, 4.5)) +
    scale_x_continuous(limits = c(0, 4.5)) +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(override.aes = list(size=3)))
}


plot_distinct_gene <- function(data, cond1, cond2, gene_expr_level) {
  label1 = paste0('<span style="color:dodgerblue4">', str_to_title(gsub("_", "\n", cond1)), '</span>')
  label2 = paste0('<span style="color:darkorange4">', str_to_title(gsub("_", "\n", cond2)), '</span>')
  facet_labels = c(label1, label2, "Up-regulated", "Down-regulated")
  names(facet_labels) = c(cond1, cond2, "UP" ,"DOWN")
  
  data %>% 
    dplyr::filter(gene_sig %in% c("UP", "DOWN")) %>% 
    dplyr::filter(enriched_in %in% c(cond1, cond2)) %>% 
    dplyr::select(-all_of(cond1), -all_of(cond2)) %>% 
    dplyr::filter(baseMean > gene_expr_level) %>% 
    dplyr::mutate(label = ifelse(str_to_lower(TF) != str_to_lower(SYMBOL), paste0(TF, "(", SYMBOL, ")"), TF),
                  label = fct_reorder(label, log2FoldChange)) %>% 
    ggplot(aes(y = log2FoldChange, x = label )) +
    geom_col(size = 0.5, position = position_dodge(), fill = "black") +
    theme_classic() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 5))+
    facet_wrap(~enriched_in, scales = "free", labeller = as_labeller(facet_labels)) +
    coord_flip() +
    theme(legend.position = "bottom", strip.text.x = element_markdown() , axis.text.y = element_text(size = 8)) +
    labs(x = "Enriched TFs" , y = "Gene Expression log2(P60/P7)") 
  
  
}

plot_pathways <-function(data, cond1, cond2, gene_expr_level){
  label1 = paste0('<span style="color:dodgerblue4">', str_to_title(gsub("_", " ", cond1)), '</span>')
  label2 = paste0('<span style="color:darkorange4">', str_to_title(gsub("_", " ", cond2)), '</span>')
  
  # pull out significant and distint TFs
  dat = data %>% 
    dplyr::filter(gene_sig %in% c("UP", "DOWN")) %>% 
    dplyr::filter(enriched_in %in% c(cond1, cond2)) %>% 
    dplyr::filter(baseMean > gene_expr_level) 
  
  # pull out gene lists
  genes_1 = mapIds(org.Mm.eg.db, keys =dat$SYMBOL[dat$enriched_in == cond1] , column ="ENTREZID", keytype="SYMBOL")
  genes_2 = mapIds(org.Mm.eg.db, keys =dat$SYMBOL[dat$enriched_in == cond2] , column ="ENTREZID", keytype="SYMBOL")
  
  # erichment analysis
  kk2 <- enrichKEGG(gene = genes_2, organism = "mmu", pvalueCutoff = 0.05)
  kk1 <- enrichKEGG(gene = genes_1, organism = "mmu", pvalueCutoff = 0.05)
  
  # plotting
  facet_labels = c(label1, label2)
  names(facet_labels) = c(cond1, cond2)
  
  bind_rows(as.data.frame(kk1) %>% dplyr::mutate(id = cond1),
            as.data.frame(kk2) %>% dplyr::mutate(id = cond2)) %>% 
    dplyr::mutate(GeneRatio_frac = as.numeric(str_extract(GeneRatio, "[^/]+")) / as.numeric(str_extract(GeneRatio, "(?<=/).+"))) %>% 
    dplyr::mutate(Description = fct_reorder(Description, GeneRatio_frac),
                  GeneRatio = fct_reorder(GeneRatio, GeneRatio_frac),
                  Count = round(Count)) %>% 
    ggplot(aes(x = GeneRatio, y = Description, size = Count, color = p.adjust))+
    geom_point() +
    facet_grid(id~., scales = "free_y",space = "free_y", labeller = as_labeller(facet_labels))+
    theme_classic() +
    theme(plot.title = element_markdown(), 
          strip.text.y = element_markdown())+
    scale_y_discrete(labels = function(x) str_wrap(x, width = 20))+
    scale_color_viridis_c(end = .9)+
    labs(y = "")
  
  
}


# >panel a ----------------------------------------------------------------
a = rasterGrob(readPNG("../figures/bart_workflow.png"), interpolate = T)

# > panel b ---------------------------------------------------------------
b = plot_distinct(early_actVrep, "early_activating", "early_repressive", "Distinct \nActivating", "Distinct \nRepressive" )

# > panel c ---------------------------------------------------------------
c = plot_distinct(late_actvrep, "late_activating", "late_repressive", "Distinct \nActivating", "Distinct \nRepressive" )

# > panel d ---------------------------------------------------------------

d = plot_distinct(actvrep, "activating", "repressive", "Distinct \nActivating", "Distinct \nRepressive" )
 

# > panel e ---------------------------------------------------------------
gene_expr_filt = 100

e1 = plot_distinct(earlyvlate, "early", "late", "Distinct Early", "Distinct Late" )
e2 = plot_distinct_gene(earlyvlate, "early", "late", gene_expr_filt )
e3 = plot_pathways(earlyvlate, "early", "late", gene_expr_filt )


#  > fig ------------------------------------------------------------------

png("../figures/poster_fig_6.png", height = 8, width = 10, res = 300, units = "in")
plot_grid(a,
          plot_grid(b,c,d, nrow = 1, labels =LETTERS[2:4]),
          plot_grid(e1,e2,e3 , nrow = 1, labels = c("Ei", "Eii", "Eiii"), rel_widths = c(1,1,1.5)),
          nrow = 3,
          labels = c("A"),
          scale = c(.9, 1, 1),
          rel_heights = c(.25, .375, .375))
dev.off()


