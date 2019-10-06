## Average methylation level over all CpGs in each sample	 
#TE_landscape/WGBS/all_CpG_Meth_average_new.txt	
python ../bin/TE_landscape/average_methylation_overall.py all_CpG_Meth.bed WGBS_samples.txt all_CpG_Meth_distribution_new.txt

## Average methylation level of CpG overlapping TEs	
#TE_landscape/WGBS/CpG_TE_Meth_average.txt	

## Average methylation level of CpG not overlapping TEs	
#TE_landscape/WGBS/CpG_noTE_Meth.bed
bedtools intersect -a all_CpG_Meth.bed -b rmsk_TEother_merge.txt -v > CpG_noTE_Meth.bed &

#TE_landscape/WGBS/CpG_noTE_Meth_average.txt	

## Average methylation level of CpGs overlapping merged Refseq features
## Split by feature to reduce memory
awk '{print $0 > $1}' WGBS/feature_CpG_Meth.bed
for file in refseq_*; do echo $file; python ~/bin/TE_landscape/average_methylation_overall.py $file sample_lists/WGBS_samples.txt $file\_average.txt feature; done 
python ~/bin/TE_landscape/average_methylation_overall.py intergenic sample_lists/WGBS_samples.txt refseq_intergenic_average.txt feature
for file in refseq_*_average.txt; do feature=$( basename "$file" _average.txt ); awk -v OFS='\t' -v feature=$feature '{print feature, $0}' $file >> WGBS/feature_CpG_Meth_average.txt; done
