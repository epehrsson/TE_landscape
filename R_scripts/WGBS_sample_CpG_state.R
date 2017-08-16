# Load number of CpGs in each methylation state by sample
# See 7/24/2017, 8/1/2017, 8/7/2017

# All CpGs
all_CpG_meth = read.table("WGBS/all_CpG_Meth_states.txt",sep=" ")
rownames(all_CpG_meth) = metadata[which(!is.na(metadata$WGBS)),]$Sample
all_CpG_meth = all_CpG_meth[,2:5]
colnames(all_CpG_meth) = meth_states

# CpGs in TEs
TE_CpG_meth = read.table("WGBS/CpG_TE_Meth_states.txt",sep=" ")
rownames(TE_CpG_meth) = metadata[which(!is.na(metadata$WGBS)),]$Sample
TE_CpG_meth = TE_CpG_meth[,2:5]
colnames(TE_CpG_meth) = meth_states
TE_CpG_meth[is.na(TE_CpG_meth)] = 0

# CpGs in RefSeq features
feature_CpG_meth = read.table("WGBS/Refseq_features/CpG_feature_Meth_states.txt",sep="\t")[,2:6]
feature_CpG_meth$Sample = rep(metadata[which(!is.na(metadata$WGBS)),]$Sample,18)
colnames(feature_CpG_meth)[1:5] = c(meth_states,"Cohort")
feature_CpG_meth[is.na(feature_CpG_meth)] = 0
feature_CpG_meth$Cohort = gsub("CpG_refseq_", "", feature_CpG_meth$Cohort)
feature_CpG_meth$Cohort = gsub("_merge_noTE_Meth", "", feature_CpG_meth$Cohort)
feature_CpG_meth$Cohort = gsub("_noTE_Meth", "", feature_CpG_meth$Cohort)