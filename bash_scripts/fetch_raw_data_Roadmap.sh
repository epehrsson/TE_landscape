# Get Roadmap raw data
# 4/18/2016, 5/3/2016, 5/5/2016, 5/19/2016, 8/25/2016, 9/12/2016, 10/27/2016, 11/3/2016, 11/4/2016, 11/7/2016, 1/6/2017, 1/18/2017, 5/10/2017, 5/25/2017, 5/29/2017, 5/30/2017, 6/5/2017, 6/12/2017, 6/21/2017, 8/2/2017

# chromHMM files
#TE_landscape/raw_data/chromHMM/all.mnemonics.bedFiles.tgz
#TE_landscape/raw_data/chromHMM/E#_15_coreMarks_mnemonics.bed [127 files]
 http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz

# WGBS
# By-chromosome WGBS fractional methylation matrices of CpG x sample 	
#/bar/mchoudhary/2chromTE/Meth/FractionalMethylation_Removed_E027_E064_Fixed_E012/

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
# hg19 RefSeq features from UCSC Table Browser
#genic_features/refseq_3UTR.txt
#genic_features/refseq_5UTR.txt
#genic_features/refseq_coding_exon.txt
#genic_features/refseq_exons.txt
#genic_features/refseq_introns.txt
#genic_features/refseq_genes.txt

# Regions 2000bp upstream of RefSeq gene TSS
#genic_features/refseq_up2000.txt

# hg19 CpG islands from UCSC Table Browser
#genic_features/cpgIslandExtUnmasked.txt

# GENCODE 
# GENCODE v19 comprehensive genes list, 2000bp upstream of TSS
genic_features/.gencodev19_2000up.bed.swp	
genic_features/Gencodev19_up2000.txt	

# GENCODE v19 comprehensive genes list
genic_features/Gencode_v19_genes.txt	

# CGate TE coordinates	
#TE_landscape/features/C_gate.txt	
https://sites.google.com/site/tecatalog/welcome/home

# VISTA enhancers	
#TE_landscape/features/vista_enhancers/vista_enhancers.txt	
https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?form=search&show=1&search.form=no&search.result=yes

# Segwey promoters
# hg19 promoters (Homo sapiens Regulatory Features (GRCh37.p13), chromosomes 1-22, X, Y, Feature Type: Promoter)
# Downloaded from ENSEMBL Biomart
#TE_landscape/features/Segway_promoters/hg19_promoters.txt

# Genome sizes
# File of chromosome sizes from UCSC	 
#TE_landscape/raw_data/hg19.genome	
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo"  > hg19.genome
genic_features/hg19.genome		

# Mouse chromosome sizes	
#TE_landscape/raw_data/mouse/mm9.genome	
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm9.chromInfo"  > mm9.genome

# Mappability
# Rsmk 36bp mappability file (Sundaram et al), standard chromosomes
#TE_landscape/mappability/rmsk.txt.mapability.36mer -> /bar/genomes/hg19/rmsk/rmsk.txt.mapability.36mer

# 36bp mappability, genome-wide	 
#TE_landscape/mappability/wgEncodeCrgMapabilityAlign36mer.bedGraph	
wget hgdownload.cse.ucsc.edu/gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign36mer.bw

# RNA-seq
# Normalization factors from ENCODE
#TE_landscape/raw_data/RNAseq/all.EGID.N.readlength

# TEs
# hg19 rmsk file
#TE_landscape/raw_data/rmsk.txt.gz -> /bar/genomes/hg19/rmsk/rmsk.txt.gz	

# RepeatMasker file for hg19	 
#TE_landscape/features/rmsk.txt	
/bar/genomes/hg19/rmsk/rmsk.txt.gz

# RepeatMasker file for mm9	 
#TE_landscape/raw_data/mouse/mm9_rmsk.txt -> /bar/genomes/mm9/rmsk/rmsk.txt	
/bar/genomes/mm9/rmsk/rmsk.txt

# RepeatMasker file for mm10	 
#TE_landscape/features/mouse/rmsk_mm10.txt	
cp /bar/genomes/mm10/rmsk/rmsk.txt.gz rmsk_mm10.txt.gz