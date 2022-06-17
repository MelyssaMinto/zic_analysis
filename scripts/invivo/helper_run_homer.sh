# Getting enriched motifs in Zic Data
# Author: Melyssa Minto

if [ -z "$1" ]; then 
        echo usage: $0 'Sample name missing'
    exit
fi


#export PATH=$HOME/meme/bin:$PATH  

# set up path to data & geneome 
SAMPLE=$1
DataDir=$2
GenomeFa=../../../../../genomeData/mm10/gencode/GRCm38.p6.genome.fa
SIZE=$3
FINDMOTIFS=~/GenomeTools/homer/bin/findMotifsGenome.pl
echo ...............$SAMPLE....................


# go to respective sample directory
echo setting up directory...
cd ../../results/invivo/homer_results/
mkdir $SAMPLE
cd $SAMPLE

# else use this for using the lengths of the peaks as is.
#echo $SAMPLE.bed > $SAMPLE.motif.bed

# get the 500bp surrounding peak
echo get the $SIZE bp surrounding peak...
SIZE_half=$(echo $(($SIZE / 2)))

cat $DataDir$SAMPLE.bed | awk -v var=$SIZE_half 'BEGIN{ OFS="\t";} { midPos=$2+$10; print $1, midPos-var , midPos+var; }' > $SAMPLE.$SIZE.bp.bed

# add other required columns 
echo formatting bed for homer...
# HOMER peak files should have at minimum 5 columns (separated by TABs, additional columns will be ignored):
#    Column1: Unique Peak ID
#    Column2: chromosome
#    Column3: starting position
#    Column5: ending position
#    Column5: Strand (+/- or 0/1, where 0="+", 1="-")


# add strand 
awk 'NR==1 {print}  NR>1 {printf("%s\t%s\n", $0, ".") }' $SAMPLE.$SIZE.bp.bed > $SAMPLE.motif.homer

# add unique peak ID
numPeaks=$(wc -l < $SAMPLE.motif.homer) # get the nummber of peaks
echo there are $numPeaks peaks

declare -A VAR # make a variable that holds the unique peak ID

for (( I=1; I <= $numPeaks; I++)); do 
	position=$(($I - 1))
	VAR[$position]="peak_${I}"
done

echo ${VAR[@]} | tr ' ' '\n' | sort > ids

paste ids $SAMPLE.motif.homer > $SAMPLE.homer

#using homer to get the de novo motifs
echo running homer..
mkdir results
perl $FINDMOTIFS $SAMPLE.homer mm10 results/invivo/ -size $SIZE -p 4 

#cleaning up files
#rm $SAMPLE.motif.homer
rm ids

echo ".................fin!......................"
