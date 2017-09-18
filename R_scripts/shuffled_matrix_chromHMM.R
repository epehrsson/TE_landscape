# Load shuffled TE potential matrices, chromHMM

# Number of samples each shuffled TE is in each chromHMM state
shuffled_chromHMM_potential = lapply(list.files(path="chromHMM/shuffled_TEs/",pattern="_potential.txt",full.names = TRUE),function(x) read.table(x,sep='\t',header=TRUE))
print("Loaded chromHMM matrices")
shuffled_chromHMM_potential = lapply(shuffled_chromHMM_potential,function(y) transform(y,States = apply(y[,8:22],1,function(x) sum(x > 0))))
print("Calculated total states")

# Maximum number of states in a single sample per TE
shuffled_max_intra = lapply(list.files(path="chromHMM/shuffled_TEs/",pattern="_max.txt",full.names = TRUE),function(x) read.table(x,sep='\t'))
shuffled_max_intra = lapply(shuffled_max_intra, setNames, nm =c("chromosome","start","stop","subfamily","class","family","strand","Max_states_intra"))
print("Loaded max matrices")
for (i in 1:10){
  shuffled_chromHMM_potential[[i]] = merge(shuffled_chromHMM_potential[[i]],shuffled_max_intra[[i]],by=c("chromosome","start","stop","subfamily","class","family","strand"),all.x=TRUE)
}
rm(shuffled_max_intra)
print("Combined chromHMM and max matrices")

save(shuffled_chromHMM_potential,file="R_datasets/shuffled_chromHMM.RData")