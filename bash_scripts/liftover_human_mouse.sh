# Human-mouse orthologs
# 9/14/2016, 9/15/2016, 9/16/2016, 2/8/2017, 2/9/2017, 5/10/2017, 5/11/2017, 5/29/2017, 6/12/2017

#To load liftOver
ml kentUCSC

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
