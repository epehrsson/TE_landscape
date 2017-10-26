# Binary matrix of TEs x sample indicating if TE is in genic or regular enhancer state
# 11/3/2016, 2/2/2017

#TE_landscape/chromHMM/enhancer_matrix.txt	
python ../bin/TE_landscape/create_state_matrix.py all_chromHMM_TE_sorted.txt rmsk_TE.txt enhancer_states.txt mnemonics.txt enhancer_matrix.txt			

#TE_landscape/chromHMM/enhancer_matrix_other.txt	
python ../bin/TE_landscape/create_state_matrix.py all_chromHMM_other_sorted.txt rmsk_other.txt enhancer_states.txt mnemonics.txt enhancer_matrix_other.txt
