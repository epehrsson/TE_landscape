# State switching
# 7/6/2016, 2/2/2017, 3/2/2017, 8/1/2017, 8/24/2017, 8/25/2017, 8/27/2017, 9/14/2017
# Updated 5/17/18 to 6/5/2018 with summit overlap

# Number of states per TE x sample
python ~/bin/TE_landscape/state_sharing_intra_lite.py chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt chromHMM/rmsk_TEother_chromHMM_summit_state_counts.txt 7 9

# State switching intra, chromHMM
#TE_landscape/chromHMM/state_switching/rmsk_TEother_chromHMM_intra.txt
python ~/bin/TE_landscape/state_sharing_intra2.py chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt chromHMM/chromHMM_states.txt chromHMM/state_switching/rmsk_TEother_chromHMM_intra.txt 7 9

# State switching inter, chromHMM
#TE_landscape/chromHMM/state_switching/rmsk_TEother_chromHMM_inter.txt
python ~/bin/TE_landscape/state_sharing_inter.py chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt chromHMM/chromHMM_states.txt chromHMM/state_switching/rmsk_TEother_chromHMM_inter.txt 7 9

# State switching inter, chromHMM, by class
#TE_landscape/chromHMM/state_switching/class/rmsk_TEother_chromHMM_inter_[class].txt [6 files]
while read line; do python ~/bin/TE_landscape/state_sharing_inter.py chromHMM/chromHMM_summit_$line\.txt chromHMM/chromHMM_states.txt chromHMM/state_switching/class/rmsk_TEother_chromHMM_inter_$line\.txt 7 9; done < features/TEs/class/TE_class.txt
python ~/bin/TE_landscape/state_sharing_inter.py chromHMM/chromHMM_summit_SVA.txt chromHMM/chromHMM_states.txt chromHMM/state_switching/class/rmsk_TEother_chromHMM_inter_SVA.txt 7 9
python ~/bin/TE_landscape/state_sharing_inter.py chromHMM/chromHMM_summit_Other.txt chromHMM/chromHMM_states.txt chromHMM/state_switching/class/rmsk_TEother_chromHMM_inter_Other.txt 7 9

# WGBS inter switching
#TE_landscape/WGBS/TE_WGBS_state_inter.txt		 
python ~/bin/TE_landscape/state_sharing_inter.py WGBS/TE_WGBS_state_sorted.txt WGBS/methylation_states.txt WGBS/TE_WGBS_state_inter.txt 7 9
 
# WGBS inter switching, by class
#TE_landscape/WGBS/class/[class]_WGBS_state_sorted.txt_inter.txt [6 files]		
for file in WGBS/class/*WGBS_state_sorted.txt; do python ~/bin/TE_landscape/state_sharing_inter.py $file WGBS/methylation_states.txt $file\_inter.txt 7 9; done &

# Number of TE x sample with CpGs in more than one state
python ~/bin/TE_landscape/shared_meth_CpG.py WGBS/TE_CpG_Meth_state.txt WGBS/TE_CpG_state_counts.txt
