#!/bin/bash
# Author: Emilano Sotelo & Melyssa Minto
# This script will use bedtools to find the overlap invitro zic peaks in the anchors of the 
# predicted chroamtin loops in adult cerebellum

cd ../results/invitro/diffbind_cutnrun_zic/

P56ANCHOR=../../../../sequencing_data/Yamada/combined_MAPS_peaks.txt
P4ANCHOR=../../../../sequencing_data/Reddy/P4_loops.bed

cat $P56ANCHOR | awk -v OFS='\t' '{ print $1, $2, $3 }' > anchors.bed
cat $P56ANCHOR | awk -v OFS='\t' '{ print $4, $5, $6 }' >> anchors.bed
cat $P4ANCHOR | awk -v OFS='\t' '{ print $1, $2, $3 }' >> anchors.bed

bedtools intersect -wa -wb -a anchors.bed -b ../cutnrun_zic_seacr_peaks_union.bed > zic_anchors.bed

