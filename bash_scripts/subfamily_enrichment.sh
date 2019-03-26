## Unique peaks overlapping each subfamily
while read line; do awk -v OFS='\t' -v sample=$line '{if($11 != "8_ZNF/Rpts") print > ""$4"_"$11"_"sample".txt"}' ~/TE_landscape/chromHMM/summit/TEs/rmsk_TEother_$line\_chromHMM.bed; while read subfamily; do while read state; do awk -v sample=$line -v subfam=$subfamily -v state=$state -v OFS='\t' '{a[$8, $9, $10]+=1}END{print state, sample, subfam, length(a)}' $subfamily\_$state\_$line\.txt >> subfamily_chromHMM_sample_summit.txt; done < chromHMM_states.txt; done < ~/TE_landscape/features/TEs/subfamily/subfamilies.txt; rm *_$line\.txt; done < mnemonics.txt
while read line; do awk -v OFS='\t' -v sample=$line '{if($11 == "8_ZNF/Rpts") print > ""$4"_8_ZNF.Rpts_"sample".txt"}' ~/TE_landscape/chromHMM/summit/TEs/rmsk_TEother_$line\_chromHMM.bed; while read subfamily; do awk -v sample=$line -v subfam=$subfamily -v OFS='\t' '{a[$8, $9, $10]+=1}END{print "8_ZNF/Rpts", sample, subfam, length(a)}' $subfamily\_8_ZNF.Rpts_$line\.txt >> subfamily_chromHMM_sample_summit_8.txt; done < ~/TE_landscape/features/TEs/subfamily/subfamilies.txt; rm *_$line\.txt; done < ~/TE_landscape/sample_lists/mnemonics.txt
cat subfamily_chromHMM_sample_summit.txt subfamily_chromHMM_sample_summit_8.txt > ~/TE_landscape/chromHMM/subfamily/subfamily_chromHMM_sample_summit.txt

## Members in state (CpGs)
#TE_landscape/WGBS/subfamily_CpG_state_members.txt
tail -n +2 TE_CpG_Meth_state.txt | awk -v OFS='\t' '{a[$4,$8]+=1; if($12 > 0) miss[$4,$8]+=1; if($9 > 0) hypo[$4,$8]+=1; if($11 > 0) hyper[$4,$8]+=1; if($10 > 0) inter[$4,$8]+=1;} END{for (i in a){split (i, sep, SUBSEP); print sep[1], sep[2], hypo[i], inter[i], hyper[i], miss[i];}}' - > subfamily_CpG_state_members.txt

## DNase merged TE subfamilies
#TE_landscape/DNase/rmsk_TEother_subfamily_DNase.txt
while read line; do bedtools intersect -wo -a ../TE_subfamilies/TEother_subfamily_merge.txt -b DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> rmsk_TEother_subfamily_DNase.txt ; done < DNase_samples.txt

## H3K27ac merged TE subfamilies
#TE_landscape/H3K27ac/rmsk_TEother_subfamily_H3K27ac.txt
while read line; do bedtools intersect -wo -a ../TE_subfamilies/TEother_subfamily_merge.txt -b H3K27ac_narrow_peaks/$line\-H3K27ac.narrowPeak | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> rmsk_TEother_subfamily_H3K27ac.txtÂ ; done < H3K27ac_samples.txt
