#! /bin/bash
# Author: Melyssa 
# Where does atoh1 binding intersect with zic binding 

# merge atoh1 rep1 and reop2 peaks

bedtools intersect -wa -wb -a ../../results/DiffExp_ZicChIP/P60vP7_DOWN.bed -b ../../results/mergedPeaks/atoh1.bed > ../../results/mergedPeaks/zic_atoh1_p7.bed
bedtools intersect -wa -wb -a ../../results/DiffExp_ZicChIP/P60vP7_UP.bed -b ../../results/mergedPeaks/atoh1.bed > ../../results/mergedPeaks/zic_atoh1_p60.bed
bedtools intersect -wa -wb -a ../../results/DiffExp_ZicChIP/P60vP7_NS.bed -b ../../results/mergedPeaks/atoh1.bed > ../../results/mergedPeaks/zic_atoh1_static.bed


./helper_meme.sh zic_atoh1_p7
./helper_meme.sh zic_atoh1_p60
./helper_meme.sh zic_atoh1_static

# ./helper_run_homer.sh zic_atoh1_p7 ../../mergedPeaks/ 500
# ./helper_run_homer.sh zic_atoh1_p60 ../../mergedPeaks/ 500
# ./helper_run_homer.sh zic_atoh1_static ../../mergedPeaks/ 500
