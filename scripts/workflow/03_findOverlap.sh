#!/bin/bash
# Author: Melyssa Minto
# This script will use bedtools to find the overlap in peaks (Zic, DNase, H3K27ac) in the anchors of the 
# predicted chroamtin loops in adult cerebellum

cd ../../results/mergedPeaks/

P56ANCHOR=../../sequencing_data/Yamada/combined_MAPS_peaks.txt
P4ANCHOR=../../sequencing_data/Reddy/P4_EP_loops.bed

cat $P56ANCHOR | awk -v OFS='\t' '{ print $1, $2, $3 }' > anchors.bed
cat $P56ANCHOR | awk -v OFS='\t' '{ print $4, $5, $6 }' >> anchors.bed

bedtools intersect -wa -wb -a anchors.bed -b DNase.bed > DNase_p56anchors.bed
bedtools intersect -wa -wb -a anchors.bed -b zic.bed > zic_p56anchors.bed
bedtools intersect -wa -wb -a anchors.bed -b H3K27ac.bed > H3K27ac_p56anchors.bed
bedtools intersect -wa -wb -a $P4ANCHOR -b DNase.bed > DNase_p4anchors.bed
bedtools intersect -wa -wb -a $P4ANCHOR -b zic.bed > zic_p4anchors.bed
bedtools intersect -wa -wb -a $P4ANCHOR -b H3K27ac.bed > H3K27ac_p4anchors.bed