#! /bin/bash
# Author: Melyssa 
# Where does atoh1 binding intersect with zic binding 

# merge atoh1 rep1 and reop2 peaks

# bedtools intersect -wa -wb -a ../../results/invivo/DiffExp_ZicChIP/P60vP7_DOWN.bed -b ../../invivo/results/mergedPeaks/atoh1.bed > ../../results/invivo/mergedPeaks/zic_atoh1_p7.bed
# bedtools intersect -wa -wb -a ../../results/invivo/DiffExp_ZicChIP/P60vP7_UP.bed -b ../../invivo/results/mergedPeaks/atoh1.bed > ../../results/invivo/mergedPeaks/zic_atoh1_p60.bed
# bedtools intersect -wa -wb -a ../../results/invivo/DiffExp_ZicChIP/P60vP7_NS.bed -b ../../invivo/results/mergedPeaks/atoh1.bed > ../../results/invivo/mergedPeaks/zic_atoh1_static.bed
# bedtools intersect -wa -wb -a ../../results/invivo/mergedPeaks/zic.bed -b ../../results/invivo/mergedPeaks/atoh1.bed > ../../results/invivo/mergedPeaks/zic_atoh1.bed


# ./helper_meme.sh zic_atoh1_p7 ../../mergedPeaks/
# ./helper_meme.sh zic_atoh1_p60 ../../mergedPeaks/
# ./helper_meme.sh zic_atoh1_static ../../mergedPeaks/
# 
# ./helper_run_homer.sh zic_atoh1_p7 ../../mergedPeaks/ 500
# ./helper_run_homer.sh zic_atoh1_p60 ../../mergedPeaks/ 500
# ./helper_run_homer.sh zic_atoh1_static ../../mergedPeaks/ 500


# find overlap
bedtools intersect -wa -wb -a ../../results/invivo/mergedPeaks/zic.bed -b ../../results/invivo/mergedPeaks/atoh1.bed > ../../results/invivo/mergedPeaks/zic_atoh1.bed
bedtools intersect -wa -wb -a ../../results/invivo/mergedPeaks/zic.bed -b ../../results/invivo/mergedPeaks/chd4P22.bed > ../../results/invivo/mergedPeaks/zic_chd4P22.bed
bedtools intersect -wa -wb -a ../../results/invivo/mergedPeaks/zic.bed -b ../../results/invivo/mergedPeaks/chd7P4.bed > ../../results/invivo/mergedPeaks/zic_chd7P4.bed
bedtools intersect -wa -wb -a ../../results/invivo/mergedPeaks/zic.bed -b ../../results/invivo/mergedPeaks/ctcfP4.bed > ../../results/invivo/mergedPeaks/zic_ctcfP4.bed
bedtools intersect -wa -wb -a ../../results/invivo/mergedPeaks/zic.bed -b ../../results/invivo/mergedPeaks/ctcfP22.bed > ../../results/invivo/mergedPeaks/zic_ctcfP22.bed
bedtools intersect -wa -wb -a ../../results/invivo/mergedPeaks/zic.bed -b ../../results/invivo/mergedPeaks/ctcfP56.bed > ../../results/invivo/mergedPeaks/zic_ctcfP56.bed
bedtools intersect -wa -wb -a ../../results/invivo/mergedPeaks/zic.bed -b ../../results/invivo/mergedPeaks/rad21P4.bed > ../../results/invivo/mergedPeaks/zic_rad21P4.bed
bedtools intersect -wa -wb -a ../../results/invivo/mergedPeaks/zic.bed -b ../../results/invivo/mergedPeaks/rad21P56.bed > ../../results/invivo/mergedPeaks/zic_rad21P56.bed
bedtools intersect -wa -wb -a ../../results/invivo/mergedPeaks/zic.bed -b ../../results/invivo/mergedPeaks/smc1p22.bed > ../../results/invivo/mergedPeaks/zic_smc1p22.bed
bedtools intersect -wa -wb -a ../../results/invivo/mergedPeaks/zic.bed -b ../../results/invivo/mergedPeaks/h2ap22.bed > ../../results/invivo/mergedPeaks/zic_h2ap22.bed
bedtools intersect -wa -wb -a ../../results/invivo/mergedPeaks/zic.bed -b ../../results/invivo/mergedPeaks/h3p22.bed > ../../results/invivo/mergedPeaks/zic_h3p22.bed
