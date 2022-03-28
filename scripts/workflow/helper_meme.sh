# !/bin/bash
# Author: Melyssa Minto

GENOME=../../../../../genomeData/mm10/gencode/GRCm38.p6.genome.fa
PREFIX=$1
DB=~/GenomeTools/meme-5.0.5/db/motif_databases/MOUSE/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme
DataDir=$2

cd ../../results/MEME
mkdir $PREFIX
cd $PREFIX


#create fasta for scaning
SIZE_half=250

cat $DataDir/$PREFIX.bed | awk -v var=$SIZE_half 'BEGIN{ OFS="\t";} { midPos=$2+$10; print $1, midPos-var , midPos+var; }' > $PREFIX.bp.bed

bedtools getfasta -fi ../$GENOME -bed $PREFIX.bp.bed -fo $PREFIX.fa
#bedtools getfasta -fi ../$GENOME -bed $DataDir/$PREFIX.bed -fo $PREFIX.fa

# Motif Scanning
echo "...motif scanning"

#> motif enrichment
meme -mod anr -maxw 15 -oc meme_results $PREFIX.fa
meme -objfun ce -minw 15 -oc meme_results_ce $PREFIX.fa

ame -oc ame_results  $PREFIX.fa $DB