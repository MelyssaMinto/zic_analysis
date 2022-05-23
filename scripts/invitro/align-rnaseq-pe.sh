#!/bin/bash
# Author: Emiliano Sotelo adapted from Melyssa Minto, VJ Ramesh Script
# Date: 4/8/2022
# This script will align RNA seq single-end fastqs

# usage: ./align.sh sample /absolute/path/to/fastq
# programs needed to run this: fastqc, STAR, UCSC genome tools aka Kent Utils, samtools, htseq, trimmomatic

# Check if the user input arguments <sample name> and <data directory>
if [ -z "$1" ]; then
 echo usage: $0 '<name> <datadir>'
 exit
elif [ -z "$2" ]; then
 echo 'no data dir'
 exit
fi

# set up all of my variables
NAME=$1 #<name>
RD=$2 #<datadir>
NAME1=${NAME}"_1"
NAME2=${NAME}"_2"
WD='/data/westlab/jes157/data_output/rnaseq/' # where the alignments will be saved
QC='/data/westlab/jes157/data_output/rnaseq/QC/'                 # where the fastQC results will be saved
TRIMMED='/data/westlab/jes157/trimmed/' # where the trimmed fastqa will be saved
trim='/data/westlab/genomeTools/Trimmomatic-0.38/trimmomatic-0.38.jar'
adapters='/data/westlab/genomeTools/Trimmomatic-0.38/adapters/'
REF='/data/westlab/reference/gencode_mouse_vM21/STAR_mm10_2.7.2b/'
GTF='/data/westlab/reference/gencode_mouse_vM21/gencode.vM21.annotation.gtf'

#begin pipeline
echo ".............$NAME.............."

echo "making sample directory..."
mkdir $WD$NAME


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

fastqc $NAME1'.fastq' --outdir=$QC
fastqc $NAME2'.fastq' --outdir=$QC

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

rm $TRIMMED$NAME1'.paired.fastq'
rm $TRIMMED$NAME2'.paired.fastq'

cd $RD
gzip -5 $NAME1'.fastq'
gzip -5 $NAME2'.fastq'

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

rm $TRIMMED$NAME1'.sorted.fastq'
rm $TRIMMED$NAME2'.sorted.fastq'

# Filter out reads that have aligned
#Flags:
#-F #: Do not output alignments with any bits set in # present in the FLAG field
#-b:Output in BAM format
#-h:Include the header in the output
#-o: ouput to file
echo getting accepted hits
samtools view -F 4 -b -h -o "$NAME"accepted_hits.bam "$NAME"Aligned.sortedByCoord.out.bam 
rm "$NAME"Aligned.sortedByCoord.out.bam 

echo building bam index
# index bamfile for bamcoverage 
samtools index -b "$NAME"accepted_hits.bam "$NAME"accepted_hits.bam.bai





##### This next section is for ATAC, ChIP, or even DNAse data ############

#Big to Bigwig
echo beginning bamcoverage normalization...

bamCoverage \
-b "$NAME"accepted_hits.bam \
-p 8  \
-o $NAME'_norm.bw' \
--normalizeUsing BPM \
--effectiveGenomeSize 2730871774 \
--ignoreForNormalization chrX 2>&1 | tee $WD$NAME'/bamcoverage.log'


echo beginning featureCounts...

echo getting counts with htseq
htseq-count --order=pos -a 30 --type=gene --format=bam --stranded=yes --idattr=gene_name "$NAME"accepted_hits.bam $GTF  > "$NAME"_gene.reads
htseq-count --order=pos -a 30 --type=exon --format=bam --stranded=yes --idattr=gene_name "$NAME"accepted_hits.bam $GTF  > "$NAME"_exon.reads



# now that we are done with everything we can zip up the fastq file to save space
rm ${TRIMMED}${NAME}*


# FIN!
echo "pipeline is done!"




