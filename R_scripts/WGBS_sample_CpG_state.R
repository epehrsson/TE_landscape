# Proportion of CpGs in each methylation state by sample
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

# Proportion
# Proportion of TE CpGs hypomethylated per sample
CpG_Meth = as.data.frame(rbind(colSums(all_CpG_meth),colSums(TE_CpG_meth)))
CpG_Meth$Cohort = c("Genome","TE")
CpG_Meth = as.data.frame(rbind(CpG_Meth,aggregate(data=feature_CpG_meth[,1:5],.~Cohort,sum)))
rownames(CpG_Meth) = CpG_Meth$Cohort
CpG_Meth = CpG_Meth[,1:4]
CpG_Meth = melt(as.matrix(CpG_Meth/rowSums(CpG_Meth)))
colnames(CpG_Meth) = c("Cohort","State","Mean")
CpG_Meth$Cohort = factor(CpG_Meth$Cohort,levels=c("TE","Genome","promoters","promoters_pc","promoters_nc","5UTR","5UTR_pc","5UTR_nc","coding_exon","coding_exon_pc","3UTR","3UTR_pc","3UTR_nc","exons","exons_pc","exons_nc","introns","introns_pc","introns_nc","intergenic"))
# mean(apply(TE_CpG_meth,1,function(x) x[1]/sum(x)))

# Contribution
# Proportion of CpGs in each state in TEs across all samples (with/without IMR90)
colSums(TE_CpG_meth)/colSums(all_CpG_meth)
colSums(TE_CpG_meth[rownames(TE_CpG_meth) != "E017",])/colSums(all_CpG_meth[rownames(all_CpG_meth) != "E017",])

# Contribution of each genic feature to CpGs in each state, across all samples
feature_CpG_meth_contribution = aggregate(data=feature_CpG_meth[,1:5],.~Cohort,sum)
rownames(feature_CpG_meth_contribution) =  feature_CpG_meth_contribution$Cohort
feature_CpG_meth_contribution =  feature_CpG_meth_contribution[,2:5]
feature_CpG_meth_contribution = t(apply(feature_CpG_meth_contribution,1,function(x) x/colSums(all_CpG_meth)))

# By-sample
# Proportion of hypomethylated CpGs in TEs by sample (see WGBS_sample_TE_state.R)
mean(TE_CpG_meth$Hypomethylated/all_CpG_meth$Hypomethylated)