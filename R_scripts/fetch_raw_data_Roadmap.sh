# Get Roadmap raw data
# 4/18/2016, 6/5/2017, 6/21/2017

# chromHMM files
#TE_landscape/raw_data/chromHMM/all.mnemonics.bedFiles.tgz
#TE_landscape/raw_data/chromHMM/E#_15_coreMarks_mnemonics.bed [127 files]
 http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz

# DNase narrow peak bedfiles
#TE_landscape/raw_data/DNase/DNase_narrow_peaks/E#-DNase.macs2.narrowPeak [53 files]
while read line; do wget http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$line\.gz; gunzip $line\.gz; done < DNase_peaks.txt

# H3K27ac narrow peak bedfiles
#TE_landscape/raw_data/H3K27ac/H3K27ac_narrow_peaks/E#-H3K27ac.narrowPeak [98 files]
while read line; do wget http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/$line\-H3K27ac.narrowPeak.gz; done < H3K27ac_samples.txt
