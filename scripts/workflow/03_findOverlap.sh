#!/bin/bash
# Author: Melyssa Minto
# This script will use bedtools to find the overlap in peaks (Zic, DNase, H3K27ac) in the anchors of the 
# predicted chroamtin loops in adult cerebellum

cd ../../results/mergedPeaks/

ANCHOR=../../sequencing_data/Yamada/combined_MAPS_peaks.txt


cat $ANCHOR | awk -v OFS='\t' '{ print $1, $2, $3 }' > anchors.bed
cat $ANCHOR | awk -v OFS='\t' '{ print $4, $5, $6 }' >> anchors.bed

bedtools intersect -wa -wb -a anchors.bed -b DNase.bed > DNase_anchors.bed
bedtools intersect -wa -wb -a anchors.bed -b zic.bed > zic_anchors.bed
bedtools intersect -wa -wb -a anchors.bed -b H3K27ac.bed > H3K27ac_anchors.bed
