# Average methylation over merged feature
# 5/25/2016, 12/14/2016, 6/15/2017, 8/7/2017

# Average methylation level over all CpGs in each sample	 
#TE_landscape/WGBS/all_CpG_Meth_average_new.txt	
python ../bin/TE_landscape/average_methylation_overall.py all_CpG_Meth.bed WGBS_samples.txt all_CpG_Meth_distribution_new.txt

#TE_landscape/WGBS/methylation_old/all_CpG_Meth_average.txt		
python ../bin/TE_landscape/average_methylation_overall.py all_CpG_Meth.bed WGBS_samples.txt
#TE_landscape/WGBS/methylation_old/error/all_CpG_Meth_average.txt		

# Average methylation level of CpG overlapping TEs	
#TE_landscape/WGBS/CpG_TE_Meth_average.txt	

# Average methylation level of CpG not overlapping TEs	
#TE_landscape/WGBS/CpG_noTE_Meth_average.txt	

# Average methylation level of CpGs overlapping merged Refseq features, no TEs	 
#TE_landscape/WGBS/Refseq_features/CpG_refseq_[feature]_merge_noTE_Meth_average.txt [6 files]	
#TE_landscape/WGBS/Refseq_features/CpG_refseq_[feature]_nc_merge_noTE_Meth_average.txt [5 files]		
#TE_landscape/WGBS/Refseq_features/CpG_refseq_[feature]_pc_merge_noTE_Meth_average.txt [6 files]		
for file in CpG_refseq_*_merge_noTE_Meth.bed; do feature=$(basename "$file" .bed); python ~/bin/TE_landscape/average_methylation_overall.py $file ../../sample_lists/WGBS_samples.txt $feature\_average.txt; done
#TE_landscape/WGBS/Refseq_features/CpG_refseq_intergenic_noTE_Meth_average.txt		 
python ~/bin/TE_landscape/average_methylation_overall.py CpG_refseq_intergenic_noTE_Meth.bed ../../sample_lists/WGBS_samples.txt CpG_refseq_intergenic_noTE_Meth_average.txt
