#!/bin/bash
# Author: Melyssa Minto
# This script will use bedtools to find the overlap in peaks (Zic, DNase, H3K27ac) in the anchors of the 
# predicted chroamtin loops in adult cerebellum

cd ../../results/mergedPeaks/

P56ANCHOR=../../sequencing_data/Yamada/combined_MAPS_peaks.txt
P4ANCHOR=../../sequencing_data/Reddy/P4_EP_loops_ext.bed

cat $P56ANCHOR | awk -v OFS='\t' '{ print $1, $2, $3 }' > anchors.bed
cat $P56ANCHOR | awk -v OFS='\t' '{ print $4, $5, $6 }' >> anchors.bed
cat $P4ANCHOR | awk -v OFS='\t' '{ print $1, $2, $3 }' >> anchors.bed

bedtools intersect -wa -wb -a anchors.bed -b DNase.bed > DNase_anchors.bed
bedtools intersect -wa -wb -a anchors.bed -b zic.bed > zic_anchors.bed
bedtools intersect -wa -wb -a anchors.bed -b H3K27ac.bed > H3K27ac_anchors.bed
