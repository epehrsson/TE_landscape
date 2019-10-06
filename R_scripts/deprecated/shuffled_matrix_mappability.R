# Add shuffled TE mappablity

load("R_datasets/shuffled_TEs_CpG.RData")

# Mappability per TE
shuffled_map = lapply(list.files(path="mappability/shuffled/",pattern="rmsk_TE_shuffle_",full.names = TRUE),function(x) read.table(x,sep='\t'))
shuffled_map = lapply(shuffled_map, setNames, nm =c("chromosome","start","stop","subfamily","class","family","strand","mappability"))
print("Loaded mappability matrices")

for (i in 1:10){
  rmsk_TE_shuffled[[i]] = merge(rmsk_TE_shuffled[[i]],shuffled_map[[i]],by=c("chromosome","start","stop","subfamily","class","family","strand"),all.x=TRUE)
}
rm(shuffled_map)
rm(i)
print("Added mappability per TE")

save(rmsk_TE_shuffled,file="R_datasets/shuffled_TEs_map.RData")