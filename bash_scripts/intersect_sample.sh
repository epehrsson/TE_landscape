#!/bin/bash

sample=$1
TEs=$2
suffix=${TEs%.txt}

#Pull out TEs in that sample
grep $sample $TEs > $suffix\_$sample\.txt
echo "Pulled out TEs for sample: "$sample

#Get 10kb bins for those TEs
bedtools intersect -wo -a $suffix\_$sample\.txt -b rmsk_TE_10kb.bed -f 1 -r -s > $suffix\_$sample\_temp.bed
awk -v OFS="\t" '{print $18,$19,$20}' $suffix\_$sample\_temp.bed > $suffix\_$sample\_10kb.bed
echo "Got 10kb regions"

#Split 10kb bins into 50bp
python /bar/epehrsson/bin/TE_landscape/split_regions_to_bins.py $suffix\_$sample\_10kb.bed $suffix\_$sample\_10kb_bins.bed 50
echo "Split 10kb regions into 50bp bins"

#Get histone files and convert to bedGraph
wget http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/$sample\-H3K4me1.fc.signal.bigwig
wget http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/$sample\-H3K4me3.fc.signal.bigwig
wget http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/$sample\-H3K27me3.fc.signal.bigwig
wget http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/$sample\-H3K36me3.fc.signal.bigwig
wget http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/$sample\-H3K9me3.fc.signal.bigwig
echo "Downloaded bigwig files"

bigWigToBedGraph $sample\-H3K27me3.fc.signal.bigwig $sample\-H3K27me3.fc.signal.bedGraph
bigWigToBedGraph $sample\-H3K36me3.fc.signal.bigwig $sample\-H3K36me3.fc.signal.bedGraph
bigWigToBedGraph $sample\-H3K4me1.fc.signal.bigwig $sample\-H3K4me1.fc.signal.bedGraph
bigWigToBedGraph $sample\-H3K4me3.fc.signal.bigwig $sample\-H3K4me3.fc.signal.bedGraph
bigWigToBedGraph $sample\-H3K9me3.fc.signal.bigwig $sample\-H3K9me3.fc.signal.bedGraph
echo "Converted to bedGraph"

rm $sample\-H3K27me3.fc.signal.bigwig
rm $sample\-H3K36me3.fc.signal.bigwig
rm $sample\-H3K4me1.fc.signal.bigwig
rm $sample\-H3K4me3.fc.signal.bigwig
rm $sample\-H3K9me3.fc.signal.bigwig

#Intersect bins with all histone modification files, then remove
bedtools intersect -wo -a $suffix\_$sample\_10kb_bins.bed -b $sample\-H3K27me3.fc.signal.bedGraph > $suffix\_$sample\-H3K27me3.fc.signal.bedGraph
bedtools intersect -wo -a $suffix\_$sample\_10kb_bins.bed -b $sample\-H3K36me3.fc.signal.bedGraph > $suffix\_$sample\-H3K36me3.fc.signal.bedGraph
bedtools intersect -wo -a $suffix\_$sample\_10kb_bins.bed -b $sample\-H3K4me1.fc.signal.bedGraph > $suffix\_$sample\-H3K4me1.fc.signal.bedGraph
bedtools intersect -wo -a $suffix\_$sample\_10kb_bins.bed -b $sample\-H3K4me3.fc.signal.bedGraph > $suffix\_$sample\-H3K4me3.fc.signal.bedGraph
bedtools intersect -wo -a $suffix\_$sample\_10kb_bins.bed -b $sample\-H3K9me3.fc.signal.bedGraph > $suffix\_$sample\-H3K9me3.fc.signal.bedGraph
echo "Intersected with histone modifications"

rm $suffix\_$sample\.txt
rm $suffix\_$sample\_temp.bed
rm $suffix\_$sample\_10kb.bed
rm $sample\-H3K27me3.fc.signal.bedGraph
rm $sample\-H3K36me3.fc.signal.bedGraph
rm $sample\-H3K4me1.fc.signal.bedGraph
rm $sample\-H3K4me3.fc.signal.bedGraph
rm $sample\-H3K9me3.fc.signal.bedGraph
rm $suffix\_$sample\_10kb_bins.bed
