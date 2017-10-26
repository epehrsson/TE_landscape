# CpGs in state by merged feature
# 6/19/2017, 6/22/2017, 7/24/2017, 7/25/2017, 8/7/2017, 8/18/2017

# Number of CpGs in each state, all CpGs, each sample
#TE_landscape/WGBS/all_CpG_Meth_states.txt		 
awk '{for (i=4;i<=NF;i++){if($i == -1) miss[i]+=1; else if ($i < 0.3) hypo[i]+=1; else if ($i > 0.7) hyper[i]+=1; else if (($i <= 0.7) && ($i >= 0.3)) inter[i]+=1;}}; END{for (i in hyper) print i" "hypo[i]" "inter[i]" "hyper[i]" "miss[i];}' all_CpG_Meth.bed > all_CpG_Meth_states.txt

#Proportion of all CpGs in each methylation state
awk '{for (i=4;i<=NF;i++){if($i == -1) sum[i]+=1;}}; END{for (i in sum) print "for column "i" is " sum[i];}' ~/TE_landscape/all_CpG_Meth.bed
awk '{for (i=4;i<=NF;i++){if(($i != -1) && ($i < 0.3)) sum[i]+=1;}}; END{for (i in sum) print "for column "i" is " sum[i];}' ~/TE_landscape/all_CpG_Meth.bed
awk '{for (i=4;i<=NF;i++){if($i > 0.7) sum[i]+=1;}}; END{for (i in sum) print "for column "i" is " sum[i];}' ~/TE_landscape/a   ll_CpG_Meth.bed
awk '{for (i=4;i<=NF;i++){if(($i <= 0.7) && ($i >= 0.3)) sum[i]+=1;}}; END{for (i in sum) print "for column "i" is " sum[i];}   ' ~/TE_landscape/all_CpG_Meth.bed

# Number of CpGs in each state, TE CpGs, each sample
#TE_landscape/WGBS/CpG_TE_Meth_states.txt		 
awk '{for (i=4;i<=NF;i++){if($i == -1) miss[i]+=1; else if ($i < 0.3) hypo[i]+=1; else if ($i > 0.7) hyper[i]+=1; else if (($i <= 0.7) && ($i >= 0.3)) inter[i]+=1;}}; END{for (i in hyper) print i" "hypo[i]" "inter[i]" "hyper[i]" "miss[i];}' CpG_TE_Meth.bed > CpG_TE_Meth_states.txt

# Proportion of TE CpGs in each methylation state
awk '{for (i=4;i<=NF;i++){if($i == -1) sum[i]+=1;}}; END{for (i in sum) print "for column "i" is " sum[i];}' ~/TE_landscape/CpG_TE_Meth.bed
awk '{for (i=4;i<=NF;i++){if(($i != -1) && ($i < 0.3)) sum[i]+=1;}}; END{for (i in sum) print "for column "i" is " sum[i];}' ~/TE_landscape/CpG_TE_Meth.bed
awk '{for (i=4;i<=NF;i++){if($i > 0.7) sum[i]+=1;}}; END{for (i in sum) print "for column "i" is " sum[i];}' ~/TE_landscape/CpG_TE_Meth.bed
awk '{for (i=4;i<=NF;i++){if(($i <= 0.7) && ($i >= 0.3)) sum[i]+=1;}}; END{for (i in sum) print "for column "i" is " sum[i];}' ~/TE_landscape/CpG_TE_Meth.bed

# Number of hypomethylated CpGs per class per sample
#TE_landscape/WGBS/class/TE_class_hypo.txt
while read line; do awk -v OFS='\t' -v class=$line '{if($5 == class) print $0}' TE_CpG_Meth_new.bed | cut -f8- - | sort | uniq | awk -v OFS='\t' -v class=$line '{for (i=4;i<=NF;i++){if(($i != -1) && ($i < 0.3)) sum[i]+=1;}}; END{for (i in sum) print i, sum[i], class;}' -; done < TEother_class.txt >>TE_class_hypo.txt
awk -v OFS='\t' '{if(($5 == "LINE?") || ($5 == "SINE?") || ($5 == "DNA?") || ($5 == "LTR?") || ($5 == "Unknown?") || ($5 == "Unknown")) print $0}' TE_CpG_Meth_new.bed | cut -f8- - | sort | uniq | awk -v OFS='\t' '{for (i=4;i<=NF;i++){if(($i != -1) && ($i < 0.3)) sum[i]+=1;}}; END{for (i in sum) print i, sum[i], "Unconfident";}' - >> TE_class_hypo.txt

# Number of CpGs in each state, TE CpGs, by sample x class
#TE_landscape/WGBS/class_CpG_Meth_states.txt		 
while read line; do awk -v OFS='\t' -v class=$line '{if($5 == class) print $0}' TE_CpG_Meth_new.bed | cut -f8- - | sort | uniq | awk -v OFS='\t' -v class=$line '{for (i=4;i<=NF;i++){if($i == -1) miss[i]+=1; else if ($i < 0.3) hypo[i]+=1; else if ($i > 0.7) hyper[i]+=1; else if (($i <= 0.7) && ($i >= 0.3)) inter[i]+=1;}}; END{for (i in hyper) print i, hypo[i], inter[i], hyper[i], miss[i], class;}' -; done < ../features/TEs/class/TEother_class.txt >> TE_class_hypo.txt
awk -v OFS='\t' '{if(($5 == "LINE?") || ($5 == "SINE?") || ($5 == "DNA?") || ($5 == "LTR?") || ($5 == "Unknown?") || ($5 == "Unknown")) print $0}' TE_CpG_Meth_new.bed | cut -f8- - | sort | uniq | awk -v OFS='\t' '{for (i=4;i<=NF;i++){if($i == -1) miss[i]+=1; else if ($i < 0.3) hypo[i]+=1; else if ($i > 0.7) hyper[i]+=1; else if (($i <= 0.7) && ($i >= 0.3)) inter[i]+=1;}}; END{for (i in hyper) print i, hypo[i], inter[i], hyper[i], miss[i], "Unconfident";}' - >> TE_class_hypo.txt (class_CpG_Meth_states.txt)
awk -v OFS='\t' '{if(($5 == "LINE?") || ($5 == "SINE?") || ($5 == "DNA?") || ($5 == "LTR?") || ($5 == "Unknown?") || ($5 == "Unknown") || ($5 == "RC")) print $0}' TE_CpG_Meth_new.bed | cut -f8- - | sort | uniq | awk -v OFS='\t' '{for (i=4;i<=NF;i++){if($i == -1) miss[i]+=1; else if ($i < 0.3) hypo[i]+=1; else if ($i > 0.7) hyper[i]+=1; else if (($i <= 0.7) && ($i >= 0.3)) inter[i]+=1;}}; END{for (i=4;i<=NF;i++) {print i, hypo[i], inter[i], hyper[i], miss[i], "Unconfident_RC";}}' - >> class_CpG_Meth_states.txt

# Number of hypomethylated CpGs per subfamily per sample
#TE_landscape/WGBS/subfamily/TE_subfamily_hypo.txt		 
while read line; do awk -v OFS='\t' -v subfam=$line '{if($4 == subfam) print $0}' TE_CpG_Meth_new.bed | cut -f8- - | sort | uniq | awk -v OFS='\t' -v subfam=$line '{for (i=4;i<=NF;i++){if(($i != -1) && ($i < 0.3)) sum[i]+=1;}}; END{for (i in sum) print i, sum[i], subfam;}' -; done < TE_subfamilies.txt >> TE_subfamily_hypo.txt

# Number of CpGs in each state, TE CpGs, by sample x subfamily
#TE_landscape/WGBS/subfamily_CpG_Meth_states.txt	
while read line; do awk -v OFS='\t' -v subfam=$line '{if($4 == subfam) print $0}' TE_CpG_Meth_new.bed | cut -f8- - | sort | uniq | awk -v OFS='\t' -v subfam=$line '{for (i=4;i<=40;i++){if($i == -1) miss[i]+=1; else if ($i < 0.3) hypo[i]+=1; else if ($i > 0.7) hyper[i]+=1; else if (($i <= 0.7) && ($i >= 0.3)) inter[i]+=1;}}; END{for (i=4;i<=40;i++){print i, hypo[i], inter[i], hyper[i], miss[i], subfam;}}' -; done < ../features/TEs/subfamily/subfamilies.txt >> subfamily_CpG_Meth_states.txt

# Number of CpGs in each state, CpGs overlapping Refseq features, each sample
#TE_landscape/WGBS/Refseq_features/CpG_feature_Meth_states.txt
 for file in CpG_refseq_*_noTE_Meth.bed; do awk -v OFS='\t' -v feature=$(basename "$file" .bed) '{for (i=4;i<=NF;i++){if($i == -1) miss[i]+=1; else if ($i < 0.3) hypo[i]+=1; else if ($i > 0.7) hyper[i]+=1; else if (($i <= 0.7) && ($i >= 0.3)) inter[i]+=1;}}; END{for (i=4;i<=NF;i++){print i, hypo[i], inter[i], hyper[i], miss[i], feature;}}' $file >>  CpG_feature_Meth_states.txt; done
