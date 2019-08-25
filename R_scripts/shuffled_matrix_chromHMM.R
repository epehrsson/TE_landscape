# Creates matrices of the number of samples each shuffled hg19 TE (row) is annotated with each chromHMM state (column) ("shuffled_chromHMM_potential")
# For 10 iterations of shuffling

# Load matrices of the number of samples each shuffled TE is in each chromHMM state
shuffled_chromHMM_potential = lapply(seq(1,10,1),function(x) 
  read.table(paste("chromHMM/shuffled_TEs/rmsk_TE_shuffle_",x,"_chromHMM_potential.txt",sep=""),sep='\t',header=TRUE))
print("Loaded chromHMM matrices")

# Add the number of unique chromHMM states each TE is annotated with across all samples
shuffled_chromHMM_potential = lapply(shuffled_chromHMM_potential,function(y) transform(y,States = apply(y[,8:22],1,function(x) sum(x > 0))))
print("Calculated total states")

# Save list of dataframes
save(shuffled_chromHMM_potential,file="R_datasets/shuffled_chromHMM.RData")