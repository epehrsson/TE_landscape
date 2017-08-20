# Proportion of each state in TEs by sample grouping

source("R_scripts/proportion.R")

# chromHMM
# Proportion of each chromHMM state in TEs by sample
mnemonics_states_TE_proportion = melt(as.matrix(mnemonics_states_TE/mnemonics_states_genome))
colnames(mnemonics_states_TE_proportion) = c("State","Sample","Proportion")

# WGBS
# Proportion of hypomethylated CpGs in TEs by sample (see WGBS_sample_TE_state.R)
WGBS_proportion = as.data.frame(cbind(rownames(TE_CpG_meth),TE_CpG_meth$Hypomethylated/all_CpG_meth$Hypomethylated))
colnames(WGBS_proportion) = c("Sample","Proportion")
WGBS_proportion$State = rep("Hypomethylated",37)

# DNase
# Proportion of Dnase peaks overlapping TEs by sample classification
DNase_proportion = as.data.frame(cbind(as.vector(DNase_stats$Sample),DNase_stats$Total_width_in_TE/DNase_stats$Total_width))
colnames(DNase_proportion) = c("Sample","Proportion")
DNase_proportion$State = rep("DNase",53)

# H3K27ac
# Proportion of H3K27ac peaks overlapping TEs by sample classification
H3K27ac_proportion = as.data.frame(cbind(as.vector(H3K27ac_stats$Sample),H3K27ac_stats$Total_width_in_TE/H3K27ac_stats$Total_width))
colnames(H3K27ac_proportion) = c("Sample","Proportion")
H3K27ac_proportion$State = rep("H3K27ac",98)