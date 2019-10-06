# Add shuffled TE CpGs

load("R_datasets/shuffled_TEs.RData")

# CpGs per TE
shuffled_WGBS_CpG = lapply(list.files(path="WGBS/shuffled/",pattern="TE_CpG_count_",full.names = TRUE),function(x) read.table(x,sep='\t'))
shuffled_WGBS_CpG = lapply(shuffled_WGBS_CpG, setNames, nm =c(TE_coordinates[c(1:4,6,5,7)],"CpGs"))
print("Loaded CpG matrices")

for (i in 1:10){
  rmsk_TE_shuffled[[i]] = merge(rmsk_TE_shuffled[[i]],shuffled_WGBS_CpG[[i]],by=TE_coordinates[c(1:4,6,5,7)],all.x=TRUE)
  rmsk_TE_shuffled[[i]][which(is.na(rmsk_TE_shuffled[[i]]$CpGs)),]$CpGs <- 0
  rmsk_TE_shuffled[[i]]$CpGs = rmsk_TE_shuffled[[i]]$CpGs/2
}
rm(shuffled_WGBS_CpG)
rm(i)
print("Added CpGs per TE")

save(rmsk_TE_shuffled,file="R_datasets/shuffled_TEs_CpG.RData")