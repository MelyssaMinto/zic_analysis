#!/bin/bash

#Author: Melyssa Minto

# This file merges all of the Zic ChIP sites and filters out blacklist. 

cd ../../results/DiffExp_ZicChIP/

## cat all files into one
echo concatenating...
cat  ../../sequencing_data/preprocessed_data/zic_chip/*/*.narrowPeak > zic.cat

##-- sorting bed file
echo sorting...
bedtools sort -i zic.cat > zic.sorted

##-- merge file 
echo merging...
bedtools merge -i zic.sorted > zic.merged

##--remove black listed regions
echo filtering...
bedtools subtract -a zic.merged -b ../../../../genomeData/mm10/bed/mm10-blacklist.bed > zic.bed

#-- cleaning up files 
rm *.cat *.merged *.sorted
