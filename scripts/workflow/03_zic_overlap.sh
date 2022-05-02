#!/bin/bash

# This script will get the overlap of k27ac and DNase with the Zic peaks 

cd ../../results/mergedPeaks

bedtools intersect -wa -wb -a zic.bed -b DNase.bed > DNase_zic_overlap.bed
bedtools intersect -wa -wb -a zic.bed -b H3K27ac.bed > H3K27ac_zic_overlap.bed

