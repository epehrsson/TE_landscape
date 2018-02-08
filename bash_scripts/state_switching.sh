# State switching
# 7/6/2016, 2/2/2017, 3/2/2017, 8/1/2017, 8/24/2017, 8/25/2017, 8/27/2017, 9/14/2017

#TE_landscape/chromHMM/state_switching/all_chromHMM_TE_state_inter.txt		
python ../bin/TE_landscape/state_sharing_inter.py all_chromHMM_TE_sorted.txt chromHMM_states.txt all_chromHMM_TE_state_inter.txt &
#TE_landscape/chromHMM/state_switching/all_chromHMM_TE_state_intra.txt		
python ../bin/TE_landscape/state_sharing_intra2.py all_chromHMM_TE_sorted.txt chromHMM_states.txt all_chromHMM_TE_state_intra.txt &
#TE_landscape/chromHMM/state_switching/all_chromHMM_other_state_inter.txt		 
python ../bin/TE_landscape/state_sharing_inter.py all_chromHMM_other_sorted.txt chromHMM_states.txt all_chromHMM_other_state_inter.txt
#TE_landscape/chromHMM/state_switching/all_chromHMM_other_state_intra.txt		 
python ../bin/TE_landscape/state_sharing_intra2.py all_chromHMM_other_sorted.txt chromHMM_states.txt all_chromHMM_other_state_intra.txt

# Number of states per TE	 
python /bar/epehrsson/bin/TE_landscape/state_sharing_inter_lite.py all_chromHMM_other_sorted.txt

# Maximum states per TE in any sample	 
#TE_landscape/chromHMM/all_chromHMM_TE_max.txt	
python ~/bin/TE_landscape/state_sharing_intra_max.py all_chromHMM_TE_sorted.txt all_chromHMM_TE_max.txt
#TE_landscape/chromHMM/all_chromHMM_other_max.txt		 
python ~/bin/TE_landscape/state_sharing_intra_max.py all_chromHMM_other_sorted.txt all_chromHMM_other_max.txt

# chromHMM inter switching, by class	 
#TE_landscape/chromHMM/state_switching/all_chromHMM_[class]_state_inter.txt [6 files]	
python ~/bin/TE_landscape/state_sharing_inter.py all_chromHMM_Other_sorted.txt chromHMM_states.txt state_switching/all_chromHMM_Other_state_inter.txt 9 7
python ~/bin/TE_landscape/state_sharing_inter.py all_chromHMM_SVA_sorted.txt chromHMM_states.txt state_switching/all_chromHMM_SVA_state_inter.txt 9 7
while read line; do python ~/bin/TE_landscape/state_sharing_inter.py all_chromHMM_$line\_sorted.txt chromHMM_states.txt state_switching/all_chromHMM_$line\_state_inter.txt 9 7; done < ../features/TEs/class/TE_class.txt

# WGBS inter switching
#TE_landscape/WGBS/TE_WGBS_state_inter.txt		 
python ~/bin/TE_landscape/state_sharing_inter.py TE_WGBS_state_sorted.txt methylation_states.txt TE_WGBS_state_inter.txt 7 9
 
# WGBS inter switching, by class
#TE_landscape/WGBS/class/[class]_WGBS_state_sorted.txt_inter.txt [6 files]		
for file in *WGBS_state_sorted.txt; do python ~/bin/TE_landscape/state_sharing_inter.py $file ../methylation_states.txt $file\_inter.txt 7 9; done

# Number of states per TE x sample (11/30/17)
 python ~/bin/TE_landscape/state_sharing_intra_lite.py chromHMM/all_chromHMM_other_sorted.txt chromHMM/chromHMM_states.txt chromHMM/all_chromHMM_other_state_counts.txt
 python ~/bin/TE_landscape/state_sharing_intra_lite.py chromHMM/all_chromHMM_TE_sorted.txt chromHMM/chromHMM_states.txt chromHMM/all_chromHMM_TE_state_counts.txt

# Number of TE x sample with CpGs in more than one state (11/30/17)
  python ~/bin/TE_landscape/shared_meth_CpG.py WGBS/TE_CpG_Meth_state.txt WGBS/TE_CpG_state_counts.txt

