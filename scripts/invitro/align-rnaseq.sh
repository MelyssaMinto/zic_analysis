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
WD='/data/westlab/jes157/data_output/rnaseq/' # where the alignments will be saved
QC='/data/westlab/jes157/data_output/rnaseq/QC/'                 # where the fastQC results will be saved
TRIMMED='/data/westlab/jes157/trimmed/' # where the trimmed fastqa will be saved
trim='/data/westlab/genomeTools/Trimmomatic-0.38/trimmomatic-0.38.jar'
adapters='/data/westlab/genomeTools/Trimmomatic-0.38/adapters/'
REF='/data/westlab/reference/gencode_mouse_vM21/STAR_mm10_2.7.2b/'
GTF='/data/westlab/reference/gencode_mouse_vM21/gencode.vM21.annotation.gtf'

###1. start aligment
echo .............$NAME..............

## get fastqc
cd $RD
# unzip file
echo unzipping fastq
FILE=$NAME'.fastq.gz'
if [ -e $NAME'.fastq.gz' ]
then
 gunzip $NAME'.fastq.gz'
fi


##2. FastQCc
echo fastQC...
fastqc $NAME'.fastq' --outdir=$QC

##3. Trim adaptors and output to a specified folder. The name of the trimmed file is $NAME + _trimmed.fastq.gz
cd $TRIMMED
echo trimming adaptors...
java -jar $trim SE -threads 8 -phred33 $RD'/'$NAME'.fastq' $NAME'.trimmed.fastq' ILLUMINACLIP:${adapters}TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


##3b. FastQCc
echo fastQC...
fastqc $TRIMMED$NAME'.trimmed.fastq' --outdir=$QC


## Run alignment
cd $WD     # changing directory to where the alignments should be saved
mkdir $NAME # making a folder to hold the samples alignment files
cd $NAME        # changing directory into the folder created above
echo run alignment...
STAR --genomeDir $REF --runThreadN 8 --readFilesIn $TRIMMED$NAME'.trimmed.fastq' --outFileNamePrefix $NAME --limitBAMsortRAM 20000000000 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard


# Filter out reads that have aligned
#Flags:
#-F #: Do not output alignments with any bits set in # present in the FLAG field
#-b:Output in BAM format
#-h:Include the header in the output
#-o: ouput to file
echo getting accepted hits
samtools view -F 4 -b -h -o "$NAME"accepted_hits.bam "$NAME"Aligned.sortedByCoord.out.bam 

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
htseq-count -a 30 --type=gene --format=bam --stranded=yes --idattr=gene_name "$NAME"accepted_hits.bam $GTF  > "$NAME"_gene.reads
htseq-count -a 30 --type=exon --format=bam --stranded=yes --idattr=gene_name "$NAME"accepted_hits.bam $GTF  > "$NAME"_exon.reads

cd $RD
gzip -5 $NAME'.fastq'

cd $TRIMMED
rm $NAME'.trimmed.fastq'

echo END



