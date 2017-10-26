# Process WGBS
# 5/5/2016, 6/15/2017

# Bedfile of CpG methylation level across all chromsomes	 
#TE_landscape/WGBS/all_CpG_Meth.bed	
for file in /bar/mchoudhary/2chromTE/Meth/FractionalMethylation_Removed_E027_E064_Fixed_E012/chr*.fm; do output=$(basename $file); chr=${output%.fm}; awk 'NR%2==0' $file | awk -v x=$chr 'BEGIN{OFS="\t";}{print x, $1-1, $0}' - >> all_CpG_Meth.bedGraph; done
sort -k1,1V -k2,2n -k3,3n all_CpG_Meth.bedGraph > all_CpG_Meth.bed #For sorting by chromosome number - does not put chr10 after chr1

# Repeated with both strands 	 
#TE_landscape/WGBS/all_CpG_Meth.bedGraph	
for file in /bar/mchoudhary/2chromTE/Meth/FractionalMethylation_Removed_E027_E064_Fixed_E012/chr*.fm; do awk -v chr=$(basename "$file" .fm) -v OFS='\t' '{print chr, $1, $1+1, $0}' $file >> all_CpG_Meth.bedGraph; done &
cut -f1-3,5- all_CpG_Meth.bedGraph | sort -k1,1V -k2,2n -k3,3n - > all_CpG_Meth.bed