# Proportion of feature in each state
# See 4/18/2016, 4/20/2016, 4/25/2016, 4/26/2016, 4/27/2016, 5/10/2016, 5/25/2016, 6/27/2016, 9/8/2016, 9/25/2016, 9/27/2016, 2/3/2017, 2/10/2017, 3/8/2017, 5/8/2017, 5/24/2017, 5/30/2017, 6/5/2017, 7/4/2017, 7/24/2017, 8/1/2017, 8/3/2017, 8/7/2017
# Features include TEs (5/31/18)

# Matrices of bases/CpGs in state by sample and cohort

# chromHMM

# Genome
# Number of bases in each state in each sample overall
mnemonics_states_genome = read.table("chromHMM/genome/mnemonics_state.txt",sep='\t',header=TRUE,row.names=1)
mnemonics_states_genome = melt(as.matrix(mnemonics_states_genome)[2:16,])
colnames(mnemonics_states_genome) = c("State","Sample","Bases")
mnemonics_states_genome$Cohort = rep("Genome",dim(mnemonics_states_genome)[1])

# TEs
# Number of bases in each state in each sample in merged TEs
mnemonics_states_TE = read.table("chromHMM/mnemonics_TEother_merge_states.txt",sep='\t',header=TRUE,row.names=1)
mnemonics_states_TE[is.na(mnemonics_states_TE)] = 0
mnemonics_states_TE = melt(as.matrix(mnemonics_states_TE)[2:16,])
colnames(mnemonics_states_TE) = c("State","Sample","Bases")
mnemonics_states_TE$Cohort = rep("TE",dim(mnemonics_states_TE)[1])

# Refseq features, no TEs
# Number of bases in each state in each sample in merged genic features, no TEs
mnemonics_states_features = read.table("chromHMM/chromHMM_refseq_features.txt",sep='\t')
colnames(mnemonics_states_features) = c("Sample","Cohort","State","Bases")
mnemonics_states_features$Cohort = gsub("refseq_","",mnemonics_states_features$Cohort)
mnemonics_expand = expand.grid(Sample = metadata$Sample,State = chromHMM_states, Cohort = unique(mnemonics_states_features$Cohort))
mnemonics_states_features = merge(mnemonics_states_features,mnemonics_expand,all.y=TRUE,all.x=TRUE,by=c("Sample","State","Cohort"))
mnemonics_states_features[which(is.na(mnemonics_states_features$Bases)),]$Bases = 0
rm(mnemonics_expand)

# Combine into one dataframe
chromHMM_state_proportion = rbind(mnemonics_states_genome,mnemonics_states_TE,mnemonics_states_features)
chromHMM_state_proportion = ddply(chromHMM_state_proportion,.(Sample,Cohort),transform,Total=sum(Bases))

# WGBS

# All CpGs
all_CpG_meth = read.table("WGBS/all_CpG_Meth_states.txt",sep=" ")
rownames(all_CpG_meth) = metadata[which(!is.na(metadata$WGBS)),]$Sample
all_CpG_meth = all_CpG_meth[,2:5]/2
colnames(all_CpG_meth) = meth_states

# CpGs in TEs
TE_CpG_meth = read.table("WGBS/CpG_TE_Meth_states.txt",sep=" ")
rownames(TE_CpG_meth) = metadata[which(!is.na(metadata$WGBS)),]$Sample
TE_CpG_meth = TE_CpG_meth[,2:5]/2
colnames(TE_CpG_meth) = meth_states
TE_CpG_meth[is.na(TE_CpG_meth)] = 0

# CpGs in RefSeq features
feature_CpG_meth = read.table("WGBS/feature_CpG_Meth_states.txt",sep="\t")
colnames(feature_CpG_meth) = c("State","Cohort","Sample","CpGs")
feature_CpG_meth$Sample = mapvalues(feature_CpG_meth$Sample,from = seq(5,41,1), to = as.vector(metadata[which(!is.na(metadata$WGBS)),]$Sample))
feature_CpG_meth$Cohort = gsub("refseq_", "", feature_CpG_meth$Cohort)
CpG_expand = expand.grid(Sample = unique(feature_CpG_meth$Sample),State = meth_states, Cohort = unique(feature_CpG_meth$Cohort))
feature_CpG_meth = merge(feature_CpG_meth,CpG_expand,all.y=TRUE,all.x=TRUE,by=c("Sample","State","Cohort"))
rm(CpG_expand)
feature_CpG_meth[is.na(feature_CpG_meth)] = 0
feature_CpG_meth$CpGs = feature_CpG_meth$CpGs/2

# Combine into one dataframe
WGBS_state_proportion = rbind(melt(as.matrix(all_CpG_meth)),melt(as.matrix(TE_CpG_meth)))
colnames(WGBS_state_proportion) = c("Sample","State","CpGs")
WGBS_state_proportion$Cohort = rep(c("Genome","TE"),each=148)
WGBS_state_proportion = rbind(WGBS_state_proportion,feature_CpG_meth)
WGBS_state_proportion = ddply(WGBS_state_proportion,.(Sample,Cohort),transform,Total=sum(CpGs))

# DNase

# Number and width of DNase peaks per sample, overall and in TEs
DNase_stats = read.table("DNase/DNase_stats.txt",sep='\t')
colnames(DNase_stats) = c("Peaks","Sample","Total_width","Peaks_in_TE")
test = read.table("DNase/rmsk_TEother_merge_DNase_contribution.txt",sep='\t')
colnames(test) = c("Sample","Total_width_in_TE")
DNase_stats = merge(DNase_stats,test,by=c("Sample"))[,c(1:2,4,3,5)]
DNase_stats$Summit_in_TE = read.table("DNase/rmsk_TEother_DNase_summit_stats.txt",sep='\t')$V1
rm(test)

# Number of bases of DNase peak in Refseq genic features, by sample
DNase_features = read.table("DNase/refseq_features_DNase.txt",sep='\t')
colnames(DNase_features) = c("Sample","Cohort","Bases")
DNase_features$Cohort = gsub("refseq_", "", DNase_features$Cohort)
DNase_features$Total = apply(DNase_features,1,function(x) ifelse(x[1] %in% as.vector(metadata[which(metadata$chrY == "Yes"),]$Sample),
                                                                 feature_overlap[which(feature_overlap$Cohort == x[2]),]$Genome,
                                                                 feature_overlap[which(feature_overlap$Cohort == x[2]),]$Genome_noY))

# Combine into one dataframe
DNase_stats_long = melt(DNase_stats[,c("Sample","Total_width","Total_width_in_TE")],id.var=("Sample"))
colnames(DNase_stats_long)[2:3] = c("Cohort","Bases")
DNase_stats_long$Total = ifelse(DNase_stats_long$Cohort == "Total_width",
                                ifelse(DNase_stats_long$Sample %in% as.vector(metadata[which(metadata$chrY == "Yes"),]$Sample),GENOME_WIDTH,GENOME_WIDTH_noY),
                                ifelse(DNase_stats_long$Sample %in% as.vector(metadata[which(metadata$chrY == "Yes"),]$Sample),MERGED_TE_WIDTH,MERGED_TE_WIDTH_noY))
DNase_stats_long = rbind(DNase_stats_long,DNase_features)
DNase_stats_long$State = rep("DNase",dim(DNase_stats_long)[1])

# H3K27ac 

# Number and width of H3K27ac peaks per sample, overall and in TEs
H3K27ac_stats = read.table("H3K27ac/H3K27ac_stats.txt",sep='\t',header=TRUE)
test = read.table("H3K27ac/rmsk_TEother_merge_H3K27ac_contribution.txt",sep='\t')
colnames(test) = c("Sample","Total_width_in_TE")
H3K27ac_stats = merge(H3K27ac_stats,test,by=c("Sample"))
H3K27ac_stats$Summit_in_TE = read.table("H3K27ac/rmsk_TEother_H3K27ac_summit_stats.txt",sep='\t')$V1
rm(test)

# Number of bases of H3K27ac peak in Refseq genic features, by sample
H3K27ac_features = read.table("H3K27ac/refseq_features_H3K27ac.txt",sep='\t')
colnames(H3K27ac_features) = c("Sample","Cohort","Bases")
H3K27ac_features$Cohort = gsub("refseq_", "", H3K27ac_features$Cohort)
H3K27ac_features$Total = apply(H3K27ac_features,1,function(x) ifelse(x[1] %in% as.vector(metadata[which(metadata$chrY == "Yes"),]$Sample),
                                                                     feature_overlap[which(feature_overlap$Cohort == x[2]),]$Genome,
                                                                     feature_overlap[which(feature_overlap$Cohort == x[2]),]$Genome_noY))

# Combine into one dataframe
H3K27ac_stats_long = melt(H3K27ac_stats[,c("Sample","Total_width","Total_width_in_TE")],id.var=("Sample"))
colnames(H3K27ac_stats_long)[2:3] = c("Cohort","Bases")
H3K27ac_stats_long$Total = ifelse(H3K27ac_stats_long$Cohort == "Total_width",
                                  ifelse(H3K27ac_stats_long$Sample %in% as.vector(metadata[which(metadata$chrY == "Yes"),]$Sample),GENOME_WIDTH,GENOME_WIDTH_noY),
                                  ifelse(H3K27ac_stats_long$Sample %in% as.vector(metadata[which(metadata$chrY == "Yes"),]$Sample),MERGED_TE_WIDTH,MERGED_TE_WIDTH_noY))

H3K27ac_stats_long = rbind(H3K27ac_stats_long,H3K27ac_features)
H3K27ac_stats_long$State = rep("H3K27ac",dim(H3K27ac_stats_long)[1])

# Combine all metrics
test = WGBS_state_proportion
colnames(test)[3] = "Bases"
all_state_proportion = rbind(chromHMM_state_proportion,test,DNase_stats_long,H3K27ac_stats_long)
rm(test)

all_state_proportion$Cohort = gsub("Total_width_in_TE","TE",all_state_proportion$Cohort)
all_state_proportion$Cohort = gsub("Total_width","Genome",all_state_proportion$Cohort)
all_state_proportion = all_state_proportion[which(!(all_state_proportion$Cohort %in% c("genes","genes_nc","genes_pc"))),]
all_state_proportion$Mark = factor(c(rep("chromHMM",38100),rep("WGBS",2960),rep("DNase",1060),rep("H3K27ac",1960)),levels=c("chromHMM","WGBS","DNase","H3K27ac"))
all_state_proportion$Proportion = all_state_proportion$Bases/all_state_proportion$Total

# Add metadata
all_state_proportion = merge(all_state_proportion,metadata[,c(1,4:9)],by="Sample")

all_state_proportion = split_coding(all_state_proportion,c("TE","Genome"))

rm(list=c("mnemonics_states_features","feature_CpG_meth","DNase_features","H3K27ac_features","chromHMM_state_proportion","WGBS_state_proportion","DNase_stats_long","H3K27ac_stats_long"))


# Proportion (across samples)
## Total feature width per sample: chromHMM/WGBS: total bases/CpGs annotated per feature; DNase/H3K27ac: total feature width
combined_proportion = merge(ddply(all_state_proportion,.(Mark,State,Cohort,Feature,Coding),summarise,Bases=sum(Bases)),
                               ddply(unique(all_state_proportion[,c("Mark","Sample","Cohort","Total")]),.(Mark,Cohort),summarise,Total=sum(Total)),by=c("Mark","Cohort"),all.x=TRUE)
combined_proportion$Proportion = combined_proportion$Bases/combined_proportion$Total
combined_proportion$State = factor(combined_proportion$State,levels=names(all_state_labels)[1:21])


# Contribution (across samples)
contribution = ddply(dcast(all_state_proportion[which(all_state_proportion$Cohort %in% c("Genome","TE")),1:4],Sample+State~Cohort,value.var="Bases"),.(State),summarise,Genome=sum(Genome),TE=sum(TE))
contribution$State = factor(contribution$State,levels=c(levels(contribution$State),"Total"))
contribution[22,] = list("Total",sum(contribution[1:15,]$Genome),sum(contribution[1:15,]$TE))
contribution$Proportion = contribution$TE/contribution$Genome

## Composite contribution
contribution_composite = as.data.frame(rbind(colSums(contribution[1:8,2:3]),colSums(contribution[9:15,2:3]),
                                             colSums(contribution[c(1:3,6:7),2:3]),colSums(contribution[4:5,2:3]),
                                             colSums(contribution[c(9,13:14),2:3]),colSums(contribution[10:12,2:3])))
rownames(contribution_composite) = c("Active","Inactive","Active Regulatory","Transcribed","Repressed","Poised Regulatory")
contribution_composite$Proportion = contribution_composite$TE/contribution_composite$Genome


# By sample proportion of each state in TEs
by_sample_all = dcast(all_state_proportion[which(all_state_proportion$Cohort %in% c("Genome","TE")),1:4],Sample+State~Cohort,value.var="Bases")
by_sample_all$Proportion = by_sample_all$TE/by_sample_all$Genome
by_sample_all = merge(by_sample_all,metadata[,c(1,4:9)],by="Sample")
