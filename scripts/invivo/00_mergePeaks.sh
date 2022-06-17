#!/bin/bash

#Author: Melyssa Minto

# This file merges all of the DHS sites and filters out blacklist. 
BLACKLIST=../../../../../genomeData/mm10/bed/mm10-blacklist.bed

cd ../../results/invivo/mergedPeaks/
files=('../../sequencing_data/preprocessed_data/DNase/*{1..3}/*.narrowPeak' '../../sequencing_data/preprocessed_data/DNase/*{7..9}/*.narrowPeak' '../../sequencing_data/preprocessed_data/DNase/*/*.narrowPeak' '../../sequencing_data/preprocessed_data/H3K27ac/*{1..2}/*.narrowPeak' '../../sequencing_data/preprocessed_data/H3K27ac/*{3..4}/*.narrowPeak' '../../sequencing_data/preprocessed_data/H3K27ac/*/*.narrowPeak' '../../sequencing_data/preprocessed_data/zic_chip/*{1..2}/*.narrowPeak' '../../sequencing_data/preprocessed_data/zic_chip/*{3..4}/*.narrowPeak' '../../sequencing_data/preprocessed_data/zic_chip/*/*.narrowPeak' '../../sequencing_data/preprocessed_data/atoh1_chip/*/*.narrowPeak' '../../sequencing_data/preprocessed_data/ctcf_p56_chip/*/*.narrowPeak' '../../sequencing_data/preprocessed_data/ctcf_p22_chip/*/*.narrowPeak' '../../sequencing_data/preprocessed_data/ctcf_p4_chip/*/*.narrowPeak' '../../sequencing_data/preprocessed_data/rad21_p56_chip/*/*.narrowPeak' '../../sequencing_data/preprocessed_data/rad21_p4_chip/*/*.narrowPeak' '../../sequencing_data/preprocessed_data/chd4_p22_chip/*/*.narrowPeak' '../../sequencing_data/preprocessed_data/chd7_p4_chip/*/*.narrowPeak' '../../sequencing_data/preprocessed_data/H2A_p22_chip/*/*.narrowPeak' '../../sequencing_data/preprocessed_data/H3_p22_chip/*/*.narrowPeak' '../../sequencing_data/preprocessed_data/smc1_p22_chip/*/*.narrowPeak')


names=(DNaseP7 DNaseP60 DNase H3K27acP7 H3K27acP60 H3K27ac ZicP7 ZicP60 zic atoh1 ctcfP56 ctcfP22 ctcfP4 rad21P56 rad21P4 chd4P22 chd7P4 h2ap22 h3p22 smc1p22)

echo "Num files "
echo "${!files[@]}"
i=0

for file in "${files[@]}"; do
  ## cat all files into one
  echo concatenating "${names[i]}"...
  echo files: ${file}
  
  eval "cat $file" > "${names[i]}".cat
  

  ##-- sorting bed file
  echo sorting "${names[i]}"...
  bedtools sort -i  "${names[i]}".cat >  "${names[i]}".sorted
  
  ##-- merge file 
  echo merging "${names[i]}"...
  bedtools merge -i  "${names[i]}".sorted >  "${names[i]}".merged
  
  ##--remove black listed regions
  echo filtering "${names[i]}"...
  bedtools subtract -a  "${names[i]}".merged -b $BLACKLIST >  "${names[i]}".bed
  
  ((i++))
  

done

#-- cleaning up files
rm *.cat *.merged *.sorted

