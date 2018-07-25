# Combine multiple epigenetic marks
# 1/30/2017, 1/31/2017, 2/9/2017, 5/15/2017, 5/16/2017, 6/6/2017, 6/7/2017, 7/4/2017, 7/20/2017, 7/21/2017, 7/27/2017, 7/28/2017, 10/5/2017
# Updated 5/7/18 and 5/22/18 with pandas script

# Combined matrix with each TE x sample and chromHMM, WGBS, DNase, H3K27ac, and RNA states, by sample
ln -s /bar/epehrsson/TE_landscape/chromHMM/rmsk_TEother_chromHMM_summit_sorted.txt chromHMM/.
ln -s /bar/epehrsson/TE_landscape/WGBS/TE_WGBS_state_sorted.txt WGBS/.
ln -s /bar/epehrsson/TE_landscape/DNase/rmsk_TEother_DNase_summit.txt DNase/.
ln -s /bar/epehrsson/TE_landscape/H3K27ac/rmsk_TEother_H3K27ac_summit.txt H3K27ac/.
ln -s /bar/epehrsson/TE_landscape/RNAseq/rmsk_TE_rpkm.txt RNA/.

awk '{print>"chromHMM/"$8}' rmsk_TEother_chromHMM_summit_sorted.txt &
awk '{print>"WGBS/"$8}' TE_WGBS_state_sorted.txt &
awk '{print>"H3K27ac/"$8}' rmsk_TEother_H3K27ac_summit.txt &
awk '{print>"DNase/"$8}' rmsk_TEother_DNase_summit.txt &
awk '{print>"RNA/"$8}' rmsk_TE_rpkm.txt & #From R

python combine_marks_pandas.py ~/TE_landscape/sample_lists/mnemonics.txt &

# Count number of TEs x sample in each state combination (all, unique)
#TE_landscape/compare_marks/combine_marks_counts.txt
#TE_landscape/compare_marks/combine_marks_counts_unique.txt
aggregate_marks_pandas.py

