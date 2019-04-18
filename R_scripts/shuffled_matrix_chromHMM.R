# Load shuffled TE potential matrices, chromHMM

# Number of samples each shuffled TE is in each chromHMM state
shuffled_chromHMM_potential = lapply(seq(1,10,1),function(x) 
  read.table(paste("chromHMM/shuffled_TEs/rmsk_TE_shuffle_",x,"_chromHMM_potential.txt",sep=""),sep='\t',header=TRUE))
print("Loaded chromHMM matrices")

shuffled_chromHMM_potential = lapply(shuffled_chromHMM_potential,function(y) transform(y,States = apply(y[,8:22],1,function(x) sum(x > 0))))
print("Calculated total states")

save(shuffled_chromHMM_potential,file="R_datasets/shuffled_chromHMM.RData")