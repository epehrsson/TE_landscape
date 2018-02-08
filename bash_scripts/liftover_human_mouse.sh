# Human-mouse orthologs
# 9/14/2016, 9/15/2016, 9/16/2016, 2/8/2017, 2/9/2017, 5/10/2017, 5/11/2017, 5/29/2017, 6/12/2017

#TE_landscape/Mouse/liftover/mm9ToHg19.over.chain.gz

# hg19 to mm10 liftover chain file	 
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToMm10.over.chain.gz

# mm9

# Liftover between human and mouse TEs	 
#TE_landscape/Mouse/liftover/rmsk_TE_hg19tomm9.bed	
ml kentUCSC #To load liftOver
liftOver rmsk_TE_numbered.txt ../../genomes/hg19/chainFiles/hg19ToMm9.over.chain.gz rmsk_TE_hg19tomm9.bed out.txt -minMatch=0.1

# Human TEs that do not lift over to mouse	
#./Mouse/out.txt	

# Numbers of human TEs that lifted over to mouse	
#TE_landscape/Mouse/liftover/numbers_old.txt	
awk '{print $4}' rmsk_TE_hg19tomm9.bed > numbers.txt

# Pull out human TEs lifted over to mouse by number	
#TE_landscape/Mouse/liftover/rmsk_TE_in_mm9.txt	
python ../bin/TE_landscape/subset_file.py rmsk_TE.txt numbers.txt rmsk_TE_in_mm9.txt &

# Combine human and lifted over mouse coordinates	 
#TE_landscape/Mouse/rmsk_TE_mm9.bed	
cut -f1-3,6 rmsk_TE_hg19tomm9.bed | paste - rmsk_TE_in_mm9.txt > rmsk_TE_mm9.bed

# Lifted-over human coordinates intersected with mouse TE coordinates	 
#TE_landscape/Mouse/liftover/hg19_mm9_TE_intersect.bed	
bedtools intersect -wo -a rmsk_TE_mm9.bed -b mm9_rmsk_TE.txt > hg19_mm9_TE_intersect.bed

# Human-mouse orthologs TEs with same subfamily name	 
#TE_landscape/Mouse/liftover/hg19_mm9_TE_intersect_same.bed	
awk -v OFS="\t" '{if($8==$15)print $0}' hg19_mm9_TE_intersect.bed > hg19_mm9_TE_intersect_same.bed

# Liftover between human and mouse, other TEs	
#TE_landscape/Mouse/liftover/rmsk_other_hg19tomm9.bed	
liftOver rmsk_other_numbered.txt /bar/genomes/hg19/chainFiles/hg19ToMm9.over.chain.gz rmsk_other_hg19tomm9.bed out_other.txt -minMatch=0.1 

# Human other TEs that do not lift over to mouse	
#TE_landscape/Mouse/liftover/out_other.txt	

#TE_landscape/Mouse/liftover/numbers_other.txt	
# Numbers of human other TEs that lifted over to mouse	 
awk '{print $4}' rmsk_other_hg19tomm9.bed > numbers.txt

# Pull out human other TEs lifted over to mouse by number	 
#TE_landscape/Mouse/liftover/rmsk_other_in_mm9.txt	
python ../bin/TE_landscape/subset_file.py rmsk_other.txt numbers.txt rmsk_other_in_mm9.txt

# Combine human and lifted over mouse coordinates	 
#TE_landscape/Mouse/rmsk_other_mm9.bed	
cut -f1-3,6 rmsk_other_hg19tomm9.bed | paste - rmsk_other_in_mm9.txt > rmsk_other_mm9.bed

# Lifted-over human other TE coordinates intersected with mouse all TE coordinates	 
#TE_landscape/Mouse/liftover/hg19_mm9_other_intersect.bed	
bedtools intersect -wo -a rmsk_other_mm9.bed -b mm9_rmsk_TEother.txt > hg19_mm9_other_intersect.bed

# Human-mouse orthologs other TEs with same subfamily name	 
#TE_landscape/Mouse/liftover/hg19_mm9_other_intersect_same.bed	
awk -v OFS="\t" '{if($8==$15)print $0}' hg19_mm9_other_intersect.bed > hg19_mm9_other_intersect_same.bed

# Lifted-over human TE coordinates intersected with mouse other TE coordinates	 
#TE_landscape/Mouse/liftover/hg19_mm9_TE_other_intersect.bed	
bedtools intersect -wo -a rmsk_TE_mm9.bed -b mm9_rmsk_other.txt > hg19_mm9_TE_other_intersect.bed

# Human-mouse orthologs TEs with same subfamily name	 
#TE_landscape/Mouse/liftover/hg19_mm9_TE_other_intersect_same.bed	
awk -v OFS="\t" '{if($8==$15)print $0}' hg19_mm9_TE_other_intersect.bed > hg19_mm9_TE_other_intersect_same.bed

# All orthologous mm9-hg19 TEs	 
#TE_landscape/Mouse/liftover/hg19_mm9_TEother_intersect_same.bed	
cat hg19_mm9_*_intersect_same.bed > hg19_mm9_TEother_intersect_same.bed #288098

# mm10

# Liftover between human TEs and mouse mm10	 
#TE_landscape/Mouse/liftover/rmsk_TE_hg19tomm10.bed	
liftOver rmsk_TEother_numbered.txt ../../genomes/hg19/chainFiles/hg19ToMm10.over.chain.gz rmsk_TE_hg19tomm10.bed out_mm10.txt -minMatch=0.1

# Human TEs that do not lift over to mouse	
#TE_landscape/Mouse/liftover/out_mm10.txt	

# Numbers of human TEs that lifted over to mouse	 
#TE_landscape/Mouse/liftover/numbers.txt	
awk '{print $4}' rmsk_TE_hg19tomm10.bed > numbers.txt

# Pull out human TEs lifted over to mouse by number	 
#TE_landscape/Mouse/liftover/rmsk_TE_in_mm10.txt	
python ../bin/TE_landscape/subset_file.py rmsk_TEother.txt numbers.txt rmsk_TE_in_mm10.txt

# Combine human and lifted over mouse coordinates	 
#TE_landscape/Mouse/rmsk_TE_mm10.bed	
cut -f1-3,6 rmsk_TE_hg19tomm10.bed | paste - rmsk_TE_in_mm10.txt > rmsk_TE_mm10.bed

# Lifted-over human coordinates intersected with mouse TE coordinates	 
#TE_landscape/Mouse/liftover/hg19_mm10_TE_intersect.bed	
bedtools intersect -wo -a rmsk_TE_mm10.bed -b mm10_rmsk_TE.txt > hg19_mm10_TE_intersect.bed

# Human-mouse orthologs TEs with same subfamily name (hg19-mm10)	 
#TE_landscape/Mouse/liftover/hg19_mm10_TE_intersect_same.bed	
awk -v OFS="\t" '{if($8==$15)print $0}' hg19_mm10_TE_intersect.bed > hg19_mm10_TE_intersect_same.bed

# mm10 TEs with orthologs in humans	 
#TE_landscape/Mouse/DNase_mm10/mm10_hg19_TE_intersect_same.bed	
awk -v OFS='\t' '{print $12, $13, $14, $15, $16, $17, $18}' ~/TE_landscape/Mouse/hg19_mm10_TE_intersect_same.bed | sort | uniq > mm10_hg19_TE_intersect_same.bed

# Test to confirm I am not missing orthologs by filtering repeats (1/19/18)
 awk -v OFS='\t' '{print $6, $7, $8, $11, $12, $13, $10}' features/mouse/rmsk_mm10.txt | bedtools intersect -wo -a Mouse/rmsk_TE_mm10.bed -b - > Mouse/hg19_mm10_repeat_intersect.bed
 awk -v OFS="\t" '{if($8==$15)print $0}' Mouse/hg19_mm10_repeat_intersect.bed > Mouse/hg19_mm10_repeat_intersect_same.bed
