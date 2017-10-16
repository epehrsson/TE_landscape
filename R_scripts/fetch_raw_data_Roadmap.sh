# Get Roadmap raw data
# 4/18/2016, 5/19/2016, 10/27/2016, 11/7/2016, 1/6/2017, 5/10/2017, 6/5/2017, 6/12/2017, 6/21/2017

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

# Mouse

# chromHMM (mm9)
#TE_landscape/raw_data/mouse/chromHMM/get_chromHMM.sh
#TE_landscape/Mouse/.mouse_chromHMM.swp
#TE_landscape/raw_data/mouse/chromHMM/ENCFF#.bed [15 files]
./mouse_chromHMM/get_chromHMM.sh

# WGBS (mm10)
#TE_landscape/raw_data/mouse/WGBS/get_meth.txt
#TE_landscape/raw_data/mouse/WGBS/ENCFF#.bed [9 files]

# DNase (mm10)
#TE_landscape/raw_data/mouse/DNase/get_files.txt
#TE_landscape/raw_data/mouse/DNase/ENCFF#.bed [20 files]
xargs -n 1 curl -O -L < get_files.txt

# Other
# CGate TE coordinates	
#TE_landscape/features/C_gate.txt	
https://sites.google.com/site/tecatalog/welcome/home

# VISTA enhancers	
#TE_landscape/features/vista_enhancers/vista_enhancers.txt	
https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?form=search&show=1&search.form=no&search.result=yes
