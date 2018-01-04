# Load shuffled TE matrices

# Shuffled TEs
rmsk_TE_shuffled = lapply(list.files(path="features/shuffled_TEs/",pattern="rmsk_TE_shuffle_",full.names = TRUE),function(x) read.table(x,sep='\t'))
rmsk_TE_shuffled = lapply(rmsk_TE_shuffled, setNames, nm=c("chromosome","start","stop","subfamily","class","family","strand"))
print("Loaded shuffled matrices")

save(rmsk_TE_shuffled,file="R_datasets/shuffled_TEs.RData")