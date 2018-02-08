# Number of CpGs in feature
# 9/5/16, 8/23/2016, 8/24/2016, 2/2/2017, 2/6/2017, 6/19/2017, 6/20/2017, 7/21/2017, 8/18/2017, 8/23/2017, 8/27/2017

# Number of CpGs per class	 
#TE_landscape/WGBS/class/TE_CpG_class.txt	
awk -v OFS='\t' '{print $5, $8, $9, $10}' TE_CpG_Meth_new.bed | sort | uniq | awk -v OFS='\t' '{a[$1]+=1}END{for(i in a){print i, a[i];}}' - > TE_CpG_class.txt 
awk -v OFS='\t' '{if($5 == "LTR?" || $5 == "DNA?" || $5 == "SINE?" || $5 == "LINE?" || $5 == "Unknown?" || $5 == "Unknown") print $8, $9, $10}' TE_CpG_Meth_new.bed | sort | uniq | wc -l >> TE_CpG_class.txt
awk -v OFS='\t' '{if($5 == "LTR?" || $5 == "DNA?" || $5 == "SINE?" || $5 == "LINE?" || $5 == "Unknown?" || $5 == "Unknown" || $5 == "RC") print $8, $9, $10}' TE_CpG_Meth_new.bed | sort | uniq | wc -l >> class/TE_CpG_class.txt #Had to add name manually

#TE_landscape/WGBS/methylation_old/class_CpG_count.txt		 
awk -v OFS='\t' '{a[$5]+=$48;}END{for(i in a) {print i, a[i];}}' TE_CpG_Meth.bed > class_CpG_count.txt
#TE_landscape/WGBS/methylation_old/other_class_CpG_count.txt		 
awk -v OFS='\t' '{a[$5]+=$48;}END{for(i in a) {print i, a[i];}}' other_CpG_Meth.bed > other/other_class_CpG_count.txt

# Number of CpGs per subfamily	 
#TE_landscape/WGBS/subfamily/TE_CpG_subfamily.txt	
awk -v OFS='\t' '{print $4, $8, $9, $10}' TE_CpG_Meth_new.bed | sort | uniq | awk -v OFS='\t' '{a[$1]+=1}END{for(i in a){print i, a[i];}}' - > TE_CpG_subfamily.txt

#TE_landscape/WGBS/methylation_old/subfamily_CpG_count.txt		
awk -v OFS='\t' '{a[$4,$5,$6]+=$48;}END{for(i in a) {split(i, sep, SUBSEP); print sep[1], sep[2], sep[3], a[i];}}' TE_CpG_Meth.bed 
#TE_landscape/WGBS/methylation_old/other_subfamily_CpG_count.txt		 
awk -v OFS='\t' '{a[$4,$5,$6]+=$48;}END{for(i in a) {split(i, sep, SUBSEP); print sep[1], sep[2], sep[3], a[i];}}' other_CpG_Meth.bed > other_subfamily_CpG_count.txt

# Number of CpGs per TE	 
#TE_landscape/WGBS/TE_CpG_count.txt	
awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' TE_CpG_Meth_new.bed > TE_CpG_count.txt

# Refseq promoters
#TE_landscape/WGBS/Refseq_promoters/refseq_promoter_unique_CpG_count.txt
 awk -v OFS='\t' '{a[$1, $2, $3, $4]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], a[i];}}' refseq_promoter_unique_CpG_Meth.bed > refseq_promoter_unique_CpG_count.txt

# Shuffled TEs
#TE_landscape/WGBS/shuffled/TE_CpG_count_#.txt [10 files]
 for j in {1..10}; do awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(k in a){split(k,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[k];}}' rmsk_TE_shuffle_$j\_Meth.bed > TE_CpG_count_$j.txt; done

# Exons (1/3/18)
 awk -v OFS='\t' '{a[$1, $2, $3, $4]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], a[i];}}' WGBS/refseq_exons_unique_CpG_Meth.bed > WGBS/refseq_exons_unique_CpG_count.txt

# Mouse TEs (12/7/17)
 awk -v OFS='\t' '{a[$1, $2, $3, $4, $5, $6, $7]+=1}END{for(i in a){split(i,sep,SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], a[i];}}' Mouse/WGBS/intersect/TEs/mm10_rmsk_TE_WGBS.bed > Mouse/WGBS/intersect/TEs/mm10_rmsk_TE_CpG_count.txt

