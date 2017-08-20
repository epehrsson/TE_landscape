# Proportion of feature in each state
# See 4/18/2016, 4/20/2016, 4/25/2016, 4/26/2016, 4/27/2016, 5/10/2016, 5/25/2016, 6/27/2016, 9/8/2016, 9/25/2016, 9/27/2016, 2/3/2017, 2/10/2017, 3/8/2017, 5/8/2017, 5/24/2017, 5/30/2017, 6/5/2017, 7/4/2017, 7/24/2017, 8/1/2017, 8/3/2017, 8/7/2017

library(plyr)
library(reshape2)

# chromHMM
# Genome
# Number of bases in each state in each sample overall
mnemonics_states_genome = read.table("chromHMM/genome/mnemonics_state.txt",sep='\t',header=TRUE,row.names=1)

# Normalized by total annotated
mnemonics_states_genome_normalized = adply(mnemonics_states_genome[2:16,],1,function(x) x/mnemonics_states_genome[1,])
rownames(mnemonics_states_genome_normalized) = chromHMM_states
mnemonics_states_genome_normalized$Mean = rowMeans(mnemonics_states_genome_normalized)
mnemonics_states_genome_normalized$SD = apply(mnemonics_states_genome_normalized,1,sd)

# TEs
# Number of bases in each state in each sample in merged TEs
mnemonics_states_TE = read.table("chromHMM/mnemonics_TEother_merge_states.txt",sep='\t',header=TRUE,row.names=1)
mnemonics_states_TE[is.na(mnemonics_states_TE)] = 0

# Normalized by total annotated
mnemonics_states_TE_normalized = adply(mnemonics_states_TE[2:16,],1,function(x) x/mnemonics_states_TE[1,])
rownames(mnemonics_states_TE_normalized) = chromHMM_states
mnemonics_states_TE_normalized$Mean = rowMeans(mnemonics_states_TE_normalized)
mnemonics_states_TE_normalized$SD = apply(mnemonics_states_TE_normalized,1,sd)

# Refseq features   #UPDATE - Current intergenic d.n. include TEs OR repeats; remove genome_noTE, repeats
# Number of bases in each state in each sample in merged genic features 
#mnemonics_states_features = read.table("chromHMM/Refseq_features/chromHMM_feature_states.txt",sep='\t')
#colnames(mnemonics_states_features) = c("State","Sample","Bases","Cohort")
#mnemonics_expand = expand.grid(Sample = metadata$Sample,State = chromHMM_states, Cohort = unique(mnemonics_states_features$Cohort))
#mnemonics_states_features = merge(mnemonics_states_features,mnemonics_expand,all.y=TRUE,all.x=TRUE,by=c("Sample","State","Cohort"))
#mnemonics_states_features[which(is.na(mnemonics_states_features$Bases)),]$Bases = 0
#test = aggregate(data=mnemonics_states_features,Bases~Sample+Cohort,sum)
#mnemonics_states_features$Proportion = apply(mnemonics_states_features,1,function(x) as.numeric(x[4])/test[which(test$Sample == x[1] & test$Cohort == x[3]),]$Bases)
#rm(test)
#mnemonics_states_features_normalized = merge(aggregate(data = mnemonics_states_features,Proportion~Cohort+State,mean),aggregate(data = mnemonics_states_features,Proportion~Cohort+State,sd),by=c("Cohort","State"))
#colnames(mnemonics_states_features_normalized)[3:4] = c("Mean","SD")

# Refseq features, no TEs
# Number of bases in each state in each sample in merged genic features, no TEs
mnemonics_states_features_noTE = read.table("chromHMM/Refseq_features/chromHMM_feature_noTE_states.txt",sep='\t')
colnames(mnemonics_states_features_noTE) = c("State","Sample","Bases","Cohort")
mnemonics_expand = expand.grid(Sample = metadata$Sample,State = chromHMM_states, Cohort = unique(mnemonics_states_features_noTE$Cohort))
mnemonics_states_features_noTE = merge(mnemonics_states_features_noTE,mnemonics_expand,all.y=TRUE,all.x=TRUE,by=c("Sample","State","Cohort"))
mnemonics_states_features_noTE[which(is.na(mnemonics_states_features_noTE$Bases)),]$Bases = 0

# Normalized by total annotated
test = aggregate(data=mnemonics_states_features_noTE,Bases~Sample+Cohort,sum)
mnemonics_states_features_noTE$Proportion = apply(mnemonics_states_features_noTE,1,function(x) as.numeric(x[4])/test[which(test$Sample == x[1] & test$Cohort == x[3]),]$Bases)
rm(test)
mnemonics_states_features_noTE_normalized = merge(aggregate(data = mnemonics_states_features_noTE,Proportion~Cohort+State,mean),aggregate(data = mnemonics_states_features_noTE,Proportion~Cohort+State,sd),by=c("Cohort","State"))
colnames(mnemonics_states_features_noTE_normalized)[3:4] = c("Mean","SD")

# Combining into one dataframe
mnemonics_states_all = rbind(mnemonics_states_genome_normalized[,128:129],mnemonics_states_TE_normalized[,128:129])
mnemonics_states_all$State = rep(chromHMM_states,2)
mnemonics_states_all$Cohort = c(rep("Genome",15),rep("TE",15))
mnemonics_states_all = rbind(mnemonics_states_all,mnemonics_states_features_noTE_normalized)
mnemonics_states_all$Feature = apply(mnemonics_states_all,1,function(x) unlist(strsplit(as.character(x[4]),"_"))[1])
mnemonics_states_all$Coding = apply(mnemonics_states_all,1,function(x) unlist(strsplit(as.character(x[4]),"_"))[2])
mnemonics_states_all[which(is.na(mnemonics_states_all$Coding)),]$Coding = "all"
mnemonics_states_all$Cohort = factor(mnemonics_states_all$Cohort,levels=c("TE","Genome","genome_noTE","promoters","promoters_pc","promoters_nc","5UTR","5UTR_pc","5UTR_nc","coding_exon","coding_exon_pc","3UTR","3UTR_pc","3UTR_nc","exons","exons_pc","exons_nc","introns","introns_pc","introns_nc","intergenic"))
mnemonics_states_all$State = factor(mnemonics_states_all$State,levels=chromHMM_states[c(1:3,6:7,4:5,8,10:12,9,13:14,15)])

# WGBS proportion 
source("~/TE_landscape/R_scripts/WGBS_sample_CpG_state.R")

# Proportion of Refseq feature and TE CpGs hypomethylated per sample
CpG_Meth = as.data.frame(rbind(colSums(all_CpG_meth),colSums(TE_CpG_meth)))
CpG_Meth$Cohort = c("Genome","TE")
CpG_Meth = as.data.frame(rbind(CpG_Meth,aggregate(data=feature_CpG_meth[,1:5],.~Cohort,sum)))
rownames(CpG_Meth) = CpG_Meth$Cohort
CpG_Meth = CpG_Meth[,1:4]
CpG_Meth = melt(as.matrix(CpG_Meth/rowSums(CpG_Meth)))
colnames(CpG_Meth) = c("Cohort","State","Mean")
CpG_Meth$Cohort = factor(CpG_Meth$Cohort,levels=c("TE","Genome","promoters","promoters_pc","promoters_nc","5UTR","5UTR_pc","5UTR_nc","coding_exon","coding_exon_pc","3UTR","3UTR_pc","3UTR_nc","exons","exons_pc","exons_nc","introns","introns_pc","introns_nc","intergenic"))
# mean(apply(TE_CpG_meth,1,function(x) x[1]/sum(x)))

# DNase proportion
source("~/TE_landscape/R_scripts/DNase_overlap.R")

# Proportion of RefSeq features overlapping DNase peaks, averaged
DNase_features = merge(DNase_features,aggregate(data=mnemonics_states_features_noTE,Bases~Cohort+Sample,sum),by=c("Sample","Cohort"))
colnames(DNase_features)[3:4] = c("Total_width_in_feature","Feature_width")
DNase_features$Proportion = DNase_features$Total_width_in_feature/DNase_features$Feature_width
DNase_proportion = aggregate(data=DNase_features,Proportion~Cohort,mean)

# Proportion of TEs overlapping DNase peaks, averaged
DNase_proportion[20,] = c("Genome",mean(DNase_stats$Total_width/as.numeric(mnemonics_states_genome[1,as.vector(DNase_stats$Sample)])))
DNase_proportion[21,] = c("TE",mean(DNase_stats$Total_width_in_TE/as.numeric(mnemonics_states_TE[1,as.vector(DNase_stats$Sample)])))
DNase_proportion$Cohort = factor(DNase_proportion$Cohort,levels=c("TE","Genome","genome_noTE","promoters","promoters_pc","promoters_nc","5UTR","5UTR_pc","5UTR_nc","coding_exon","coding_exon_pc","3UTR","3UTR_pc","3UTR_nc","exons","exons_pc","exons_nc","introns","introns_pc","introns_nc","intergenic"))
DNase_proportion$Proportion = as.numeric(DNase_proportion$Proportion)

# H3K27ac proportion
source("~/TE_landscape/R_scripts/H3K27ac_overlap.R")

# Proportion of RefSeq features overlapping H3K27ac peaks, averaged
H3K27ac_features = merge(H3K27ac_features,aggregate(data=mnemonics_states_features_noTE,Bases~Cohort+Sample,sum),by=c("Sample","Cohort"))
colnames(H3K27ac_features)[3:4] = c("Total_width_in_feature","Feature_width")
H3K27ac_features$Proportion = H3K27ac_features$Total_width_in_feature/H3K27ac_features$Feature_width
H3K27ac_proportion = aggregate(data=H3K27ac_features,Proportion~Cohort,mean)

# Proportion of TEs overlapping H3K27ac peaks, averaged
H3K27ac_proportion[20,] = c("Genome",mean(H3K27ac_stats$Total_width/as.numeric(mnemonics_states_genome[1,as.vector(H3K27ac_stats$Sample)])))
H3K27ac_proportion[21,] = c("TE",mean(H3K27ac_stats$Total_width_in_TE/as.numeric(mnemonics_states_TE[1,as.vector(H3K27ac_stats$Sample)])))
H3K27ac_proportion$Cohort = factor(H3K27ac_proportion$Cohort,levels=c("TE","Genome","genome_noTE","promoters","promoters_pc","promoters_nc","5UTR","5UTR_pc","5UTR_nc","coding_exon","coding_exon_pc","3UTR","3UTR_pc","3UTR_nc","exons","exons_pc","exons_nc","introns","introns_pc","introns_nc","intergenic"))
H3K27ac_proportion$Proportion = as.numeric(H3K27ac_proportion$Proportion)