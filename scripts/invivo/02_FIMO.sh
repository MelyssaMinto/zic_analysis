# #!/bin/bash
# Author: Melyssa Minto (msm110@duke.edu)
# this script will run FIMO for motif occurances 


MOTIF=../../sequencing_data/motif_data/ZIC12_HOMOCOCO.pwm
PEAK=../../results/mergedPeaks/zic.bed
GENOME=../../../../../genomeData/mm10/gencode/GRCm38.p6.genome.fa
PREFIX=Zic



mkdir ../../results/invivo/FIMO

#create fasta for scaning
bedtools getfasta -fi $GENOME -bed $PEAK -fo '../../results/invivo/FIMO/'$PREFIX.fa

# Motif Scanning
echo "...motif scanning"
#> scan condition for target motif
fimo --o '../../results/invivo/FIMO/'$PREFIX $MOTIF '../../results/invivo/FIMO/'$PREFIX.fa


echo 'Fin!'