#!/bin/bash
# Author: Melyssa Minto (msm110@duke.edu)
# this script will run FIMO for motif occurrences of Zic1 and Zic2 motifs at ZDD genes


MOTIF=../../sequencing_data/motif_data/ZIC12_HOMOCOCO.pwm
PEAK=(../../results/invitro/beds/zdd_anchors.bed ../../results/invitro/beds/zd_anchors.bed ../../results/invitro/beds/dev_anchors.bed ../../results/invitro/beds/zdd_zicPeak.bed ../../results/invitro/beds/zd_zicPeak.bed ../../results/invitro/beds/dev_zicPeak.bed)
GENOME=../../../../../genomeData/mm10/gencode/GRCm38.p6.genome.fa
PREFIX=(zdd zd dev zdd_zic zd_zic dev_zic)



mkdir ../../results/invitro/FIMO
i=0
#create fasta for scaning
for bed in "${PEAK[@]}"; do
  name=${PREFIX[i]}

  echo  "$name"...
  echo files: ${bed}
  
  # make fastq file
  bedtools getfasta -fi $GENOME -bed $bed -fo  '../../results/invitro/FIMO/'$name.fa
  

  # Motif Scanning
  echo "...motif scanning"
  #> scan condition for target motif
  fimo --o '../../results/invitro/FIMO/'$name $MOTIF '../../results/invitro/FIMO/'$name.fa
  #> meme motif enrichment
  # ./helper_meme.sh $name'_anchors' ../../beds/
  #> homer motif enrichment
  ./helper_run_homer.sh $name'_anchors' ../../beds/ 500
  
   ((i++))
done


echo 'Fin!'