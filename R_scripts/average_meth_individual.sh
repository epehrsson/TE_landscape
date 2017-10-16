# Average methylation level over individual feature
# 5/9/2016, 12/14/16, 2/2/2017, 5/11/2017, 6/19/2017, 8/7/2017, 8/23/2017, 8/28/2017

# Average methylation level of each TE in each sample	 
#TE_landscape/WGBS/TE_CpG_Meth_new_average.txt	
python ~/bin/TE_landscape/average_methylation.py TE_CpG_Meth_new.bed ~/TE_landscape/rmsk_TEother.txt ~/TE_landscape/WGBS_samples.txt TE_CpG_Meth_new_average.txt

#TE_landscape/WGBS/methylation_old/TE_CpG_Meth_avg.bed		
python ../bin/TE_landscape/average_methylation.py TE_CpG_Meth.bed rmsk_TE.txt WGBS_samples.txt
#TE_landscape/WGBS/methylation_old/other_CpG_Meth_avg.bed		 
python ../bin/TE_landscape/average_methylation.py other_CpG_Meth.bed rmsk_other.txt WGBS_samples.txt other_CpG_Meth_avg.bed
#TE_landscape/WGBS/methylation_old/error/TE_CpG_Meth_avg.bed		

# Average methylation of unique Refseq promoters	 
#TE_landscape/WGBS/Refseq_promoters/refseq_promoter_unique_CpG_Meth_average.txt	
python ~/bin/TE_landscape/average_methylation_promoter.py refseq_promoter_unique_CpG_Meth.bed ~/genic_features/RefSeq/refseq_promoters_unique_std.txt ../../sample_lists/WGBS_samples.txt refseq_promoter_unique_CpG_Meth_average.txt

# Average methylation per shuffled TE	 
#TE_landscape/WGBS/shuffled/rmsk_TE_shuffle_#_Meth_average.txt [10 files]	
for i in {1..10}; do python ~/bin/TE_landscape/average_methylation.py rmsk_TE_shuffle_$i\_Meth.bed ../rmsk_TE_shuffle_$i\.txt ~/TE_landscape/sample_lists/WGBS_samples.txt rmsk_TE_shuffle_$i\_Meth_average.txt; done

# Mouse
# mm10 TE average methylation 
#TE_landscape/Mouse/WGBS/mm10_rmsk_TE_WGBS_avg.txt	
python ~/bin/TE_landscape/average_methylation_stranded.py mm10_rmsk_TE_WGBS.bed mm10_rmsk_TE.txt samples.txt  mm10_rmsk_TE_WGBS_avg.txt
