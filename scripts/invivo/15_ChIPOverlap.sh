#!/bin/bash
# Author: Melyssa Minto

# making chip comparison plots

DATADIR=../../sequencing_data/preprocessed_data/
P60PEAKS=../../results/invivo/DiffExp_ZicChIP/P60vP7_UP.bed
P7PEAKS=../../results/invivo/DiffExp_ZicChIP/P60vP7_DOWN.bed
STATCPEAKS=../../results/invivo/DiffExp_ZicChIP/P60vP7_STATIC.bed
ALLPEAKS=../../results/invivo/mergedPeaks/zic.bed
GENESUP=../../results/invivo/peak_gene/zicPeaks_P60genes.bed
GENESDOWN=../../results/invivo/peak_gene/zicPeaks_P7genes.bed
GENESNS=../../results/invivo/peak_gene/zicPeaks_staticgenes.bed
ACTIVATION=../../results/invivo/peak_gene/activating/activating.bed
REPRESSION=../../results/invivo/peak_gene/repressive/repressive.bed


all_chip=($DATADIR'zic_chip/SRR1557091/SRR1557091_norm.bw'  $DATADIR'zic_chip/SRR1557093/SRR1557093_norm.bw' $DATADIR'ctcf_p4_chip/SRR13371014/SRR13371014_norm.bw' $DATADIR'ctcf_p56_chip/SRR8696015/SRR8696015_norm.bw' $DATADIR'rad21_p4_chip/SRR13371020/SRR13371020_norm.bw' $DATADIR'rad21_p56_chip/SRR8696023/SRR8696023_norm.bw' $DATADIR'chd7_p4_chip/SRR13371039/SRR13371039_norm.bw' $DATADIR'chd4_p22_chip/SRR3659055/SRR3659055_norm.bw')
cohesin=($DATADIR'zic_chip/SRR1557091/SRR1557091_norm.bw'  $DATADIR'zic_chip/SRR1557093/SRR1557093_norm.bw' $DATADIR'ctcf_p4_chip/SRR13371014/SRR13371014_norm.bw' $DATADIR'ctcf_p56_chip/SRR8696015/SRR8696015_norm.bw' $DATADIR'rad21_p4_chip/SRR13371020/SRR13371020_norm.bw' $DATADIR'rad21_p56_chip/SRR8696023/SRR8696023_norm.bw' )
remodel=($DATADIR'zic_chip/SRR1557091/SRR1557091_norm.bw'  $DATADIR'zic_chip/SRR1557093/SRR1557093_norm.bw' $DATADIR'rad21_p4_chip/SRR13371020/SRR13371020_norm.bw' $DATADIR'chd7_p4_chip/SRR13371039/SRR13371039_norm.bw' $DATADIR'chd4_p22_chip/SRR3659055/SRR3659055_norm.bw')
all_chip_label=("Zic P7" "Zic P60" "CTCF P4" "CTCF P56" "Rad21 P4" "Rad21 P56" "ChD7 P4" "ChD4 P22")

echo array length is: "${#all_chip[@]}"
echo array name length is: "${#all_chip[@]}"


######### Early Zic with Early marks #############
echo all peaks...
computeMatrix scale-regions -S ${all_chip[@]} -R  $P7PEAKS $STATICPEAKS $P60PEAKS -a 1000 -b 1000 -o ../../results/invivo/ChIPProfile/all_matrix.txt.gz --skipZeros	--outFileSortedRegions ../../results/invivo/ChIPProfile/all_regions_sorted.bed || { echo failed; exit 1;}
plotHeatmap -m ../../results/invivo/ChIPProfile/all_matrix.txt.gz -out ../../figures/dev_clus_zic_hm.png --kmeans 3  --dpi 100
plotHeatmap -m ../../results/invivo/ChIPProfile/all_matrix.txt.gz -out ../../figures/dev_zic_hm.png --dpi 100


####### Peaks by genes #######
echo mapped genes
computeMatrix scale-regions -S ${all_chip[@]} -R  $GENESUP $GENESDOWN $GENESNS -a 1000 -b 1000 -o ../../results/invivo/ChIPProfile/dev_gene_matrix.txt.gz --skipZeros	--outFileSortedRegions ../../results/invivo/ChIPProfile/dev_gene_regions_sorted.bed || { echo failed; exit 1;}
plotHeatmap -m ../../results/invivo/ChIPProfile/dev_gene_matrix.txt.gz -out ../../figures/dev_gene_clus_hm.png --kmeans 3 --dpi 100
plotHeatmap -m ../../results/invivo/ChIPProfile/dev_gene_matrix.txt.gz -out ../../figures/dev_gene_hm.png --dpi 100


######## Peaks by regulation ########
echo regulation
computeMatrix scale-regions -S ${all_chip[@]} -R  $ACTIVATION $REPRESSION -a 1000 -b 1000 -o ../../results/invivo/ChIPProfile/reg_matrix.txt.gz --skipZeros	--outFileSortedRegions ../../results/invivo/ChIPProfile/reg_regions_sorted.bed ||{ echo failed; exit 1;}
plotHeatmap -m ../../results/invivo/ChIPProfile/reg_matrix.txt.gz -out ../../figures/reg_hm.png  --dpi 100
plotHeatmap -m ../../results/invivo/ChIPProfile/reg_matrix.txt.gz -out ../../figures/reg_clus_hm.png  --kmeans 3 --dpi 100

