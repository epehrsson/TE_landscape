# Add shuffled TE CpGs

load("R_datasets/shuffled_TEs.RData")

# CpGs per TE
shuffled_WGBS_CpG = lapply(list.files(path="WGBS/shuffled/",pattern="TE_CpG_count_",full.names = TRUE),function(x) read.table(x,sep='\t'))
shuffled_WGBS_CpG = lapply(shuffled_WGBS_CpG, setNames, nm =c("chromosome","start","stop","subfamily","class","family","strand","CpGs"))
print("Loaded CpG matrices")

for (i in 1:10){
  rmsk_TE_shuffled[[i]] = merge(rmsk_TE_shuffled[[i]],shuffled_WGBS_CpG[[i]],by=c("chromosome","start","stop","subfamily","class","family","strand"),all.x=TRUE)
  rmsk_TE_shuffled[[i]][which(is.na(rmsk_TE_shuffled[[i]]$CpGs)),]$CpGs <- 0
}
rm(shuffled_WGBS_CpG)
rm(i)
print("Added CpGs per TE")

save(rmsk_TE_shuffled,file="R_datasets/shuffled_TEs_CpG.RData")