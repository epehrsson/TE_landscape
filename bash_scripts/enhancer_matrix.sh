# Binary matrix of TEs x sample indicating if TE is in genic or regular enhancer state
# 11/3/2016, 2/2/2017
# Updated 6/21/18

# From files with TEs in the 6_EnhG/7_Enh state per sample (combine marks)
cat chromHMM/enhancer_matrix/* | cut -f1-8 - > chromHMM/enhancer_matrix.txt
