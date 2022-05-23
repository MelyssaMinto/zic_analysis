# !/bin/bash
# Author: Vijyendra Ramesh
# Modified: Emiliano Sotelo
# Updated: 03/09/22

#  This script will align CUT and RUN DNA 
# usage ./align_paired <read1> <read2> <data dir>

# check for correct inputs
if [ -z "$1" ]; then
  echo "usage: $0 <Read1> <Read2>  <DataDir>"
  exit
elif [ -z "$2" ]; then
  echo "havent supplied name of second read"
  exit
fi

#setting up variables
NAME1=$1 #<NAME>
NAME2=$2
RD=$3 #<DataDir>
NAME=$(echo ${NAME1%_R*}) # remove the read information from the names(i.e. sample_name_lane_R1 -> sample_name_lane)
WD='/data/westlab/jes157/data_output/cutnrun/' #<working_directory>
TRIMMED='/data/westlab/jes157/trimmed/'
trim='/data/westlab/genomeTools/Trimmomatic-0.38/trimmomatic-0.38.jar'
adapters='/data/westlab/genomeTools/Trimmomatic-0.38/adapters/'
REF='/data/westlab/reference/gencode_mouse_vM21/STAR_mm10_2.7.2b/'
GTF='/data/westlab/reference/gencode_mouse_vM21/gencode.vM21.annotation.gtf'

#begin pipeline
echo ".............$NAME.............."

echo "making sample directory..."
mkdir $WD$NAME
mkdir $WD$NAME'/QC/'

#1. unzipping files
echo "unzipping fastq file..."
cd $RD
# unziping read1
if [ -e $NAME1'.fastq.gz' ]
then
 gunzip $NAME1'.fastq.gz'
fi
# unzipping read2
if [ -e $NAME2'.fastq.gz' ]
then
 gunzip $NAME2'.fastq.gz'
fi

#2.  Get fastQC
echo "..........2. FastQC................"

fastqc $NAME1'.fastq' --outdir=$WD$NAME'/QC/'
fastqc $NAME2'.fastq' --outdir=$WD$NAME'/QC/'

#3. Trim adaptors 
echo "..........3. Trimming Adaptors......"
java -jar $trim PE \
-threads 2 \
-phred33 \
$RD'/'$NAME1'.fastq' $RD'/'$NAME2'.fastq' $TRIMMED$NAME1'.paired.fastq' $TRIMMED$NAME1'.unpaired.fastq' $TRIMMED$NAME2'.paired.fastq' $TRIMMED$NAME2'.unpaired.fastq' \
ILLUMINACLIP:$adapters'/'TruSeq3-PE-NEB.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36 2>&1 | tee $WD$NAME'/trimmonatic.log'

#4. sorting trimmed reads 
echo "...........4. Sorting Trimmed Fastqs..."
cat $TRIMMED$NAME1'.paired.fastq' | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" >$TRIMMED$NAME1'.sorted.fastq'
cat $TRIMMED$NAME2'.paired.fastq' | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" >$TRIMMED$NAME2'.sorted.fastq'

#5. Align to mm10#

echo ".............6.  Aligning Sorted Fastq to MM10 genome"

cd $WD$NAME
STAR \
--genomeDir $REF \
--runThreadN 8 \
--readFilesIn $TRIMMED$NAME1'.sorted.fastq' $TRIMMED$NAME2'.sorted.fastq' \
--outFileNamePrefix $NAME \
--limitBAMsortRAM 20000000000 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard

#7.creating normalized non-stranded bigwig with deeptools

echo ".............7.  Beginning Bamcoverage and Spike-in Normalization..."


#index bamfile for bamcoverage 
samtools index -b *Aligned.sortedByCoord.out.bam "$NAME"Aligned.sortedByCoord.out.bam.bai



bamCoverage \
-b *Aligned.sortedByCoord.out.bam \
-p 8  \
-o $NAME'_norm.bw' \
--normalizeUsing BPM \
--effectiveGenomeSize 2654621837 \
--ignoreForNormalization chrX 2>&1 | tee $WD$NAME'/bamcoverage.log'


#8.BAM to Bed Conversion
echo ".............8.  Convert sortbyCoord.bam to Bed..."

bedtools bamtobed -i "$NAME"Aligned.sortedByCoord.out.bam > sequence.bed || { echo 'ERROR: could not produce sequence.bed from .bam alignment, exiting...' ; exit 1; }

#Sorting feature file by Chromosome
echo 'Starting Filtering...'
CRMREFIX='chr'
    for (( i = 1 ; i <= 25; i++ ))
    do
    CRM=$CRMREFIX$i
    grep -w $CRM sequence.bed > "$CRM".bed  #output chromosome specific bed file
    bedtools sort -i "$CRM".bed > sorted."$CRM".bed
    done
grep -w "chrX\|chrY" sequence.bed > chrXY.bed
bedtools sort -i chrXY.bed > sorted.chrXY.bed
cat sorted.chr*bed > sorted.all.bed #combining all sorted feature files
grep -w "chrM" sequence.bed | wc -l > chrM.txt

# Filtering duplicates
echo "Filtering Duplicates"
macs2 filterdup -i sorted.all.bed --keep-dup=1 -o sequence.final.bed

rm *chr*

#9. Peak Calling Using SEACR
echo "................9.  Calling CUT and RUN  peaks with SEACR..."

module load Anaconda3/2019.10-gcb02
source activate /data/westlab/conda

chromSize='/data/westlab/reference/gencode_mouse_vM21/STAR_mm10_2.7.2b/chrNameLength.txt'
seacr='/data/westlab/genomeTools/SEACR/SEACR_1.3.sh'

bedtools genomecov -bg -scale 1.0 -i sequence.final.bed -g $chromSize > sequence.normalized.bedgraph 
bash $seacr sequence.normalized.bedgraph 0.01 non stringent "$NAME".seacr_peaks.bed


echo "................10.  Finished Peak Calling, Now Tidying Up"

#10..............Tidying Up................

cd $RD
gzip -5 $NAME1'.fastq'
gzip -5 $NAME2'.fastq'

cd $TRIMMED
rm $TRIMMED$NAME1'.paired.fastq'
rm $TRIMMED$NAME2'.paired.fastq'
rm $TRIMMED$NAME1'.unpaired.fastq'
rm $TRIMMED$NAME2'.unpaired.fastq'
rm $TRIMMED$NAME1'.sorted.fastq'
rm $TRIMMED$NAME2'.sorted.fastq'

cd $WD$NAME
rm *chr*
rm sequence*
rm sorted.all.bed

echo '............... FINISH'

