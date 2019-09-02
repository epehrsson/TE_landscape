# Creates dataframes of the length/proportion of each state within features
## mnemonics_states_genome: Total bp per chromHMM state per sample
## mnemonics_states_TE: Total bp per chromHMM state overlapping TEs per sample
## all_CpG_meth: Number of CpGs in each methylation state per sample
## TE_CpG_meth: Number of CpGs overlapping TEs in each methylation state per sample
## DNase_stats: Number and total width of DHS peaks per sample, overall and overlapping TEs
## H3K27ac_stats: Number and total width of H3K27ac peaks per sample, overall and overlapping TEs
## all_state_proportion: Length and proportion of feature in state by sample, all techniques
## combined_proportion: Proportion of feature in state, across samples
## contribution: Proportion of state within TEs, across samples
## contribution_composite: Proportion of composite chromHMM state within TEs, across samples
## by_sample_all: Proportion of state within TEs, by sample

# Load matrices of bases/CpGs in state by technique, feature, and sample

# chromHMM

# Genome
# Number of bases annotated with each state in each sample
mnemonics_states_genome = read.table("chromHMM/genome/mnemonics_state.txt",sep='\t',header=TRUE,row.names=1)
mnemonics_states_genome = melt(as.matrix(mnemonics_states_genome)[2:16,])
colnames(mnemonics_states_genome) = c("State","Sample","Bases")
mnemonics_states_genome$Cohort = rep("Genome",dim(mnemonics_states_genome)[1])

# TEs
# Number of bases annotated with each state in each sample within TEs
mnemonics_states_TE = read.table("chromHMM/mnemonics_TEother_merge_states.txt",sep='\t',header=TRUE,row.names=1)
mnemonics_states_TE[is.na(mnemonics_states_TE)] = 0
mnemonics_states_TE = melt(as.matrix(mnemonics_states_TE)[2:16,])
colnames(mnemonics_states_TE) = c("State","Sample","Bases")
mnemonics_states_TE$Cohort = rep("TE",dim(mnemonics_states_TE)[1])

# RefSeq genic features
# Number of bases annotated with each state in each sample within genic features
mnemonics_states_features = read.table("chromHMM/chromHMM_refseq_features.txt",sep='\t',col.names=c("Sample","Cohort","State","Bases"))
mnemonics_states_features$Cohort = gsub("refseq_","",mnemonics_states_features$Cohort)

# Include all sample x feature x state combinations
mnemonics_expand = expand.grid(Sample = metadata$Sample,State = chromHMM_states, Cohort = unique(mnemonics_states_features$Cohort))
mnemonics_states_features = merge(mnemonics_states_features,mnemonics_expand,all.y=TRUE,all.x=TRUE,by=c("Sample","State","Cohort"))
mnemonics_states_features[which(is.na(mnemonics_states_features$Bases)),]$Bases = 0
rm(mnemonics_expand)

# Combine into one dataframe
chromHMM_state_proportion = rbind(mnemonics_states_genome,mnemonics_states_TE,mnemonics_states_features)

# Add a column listing the total bp annotated by chromHMM per feature and sample
chromHMM_state_proportion = ddply(chromHMM_state_proportion,.(Sample,Cohort),transform,Total=sum(Bases))

# WGBS

# Genome
# Number of CpGs in each methylation state per sample
all_CpG_meth = read.table("WGBS/all_CpG_Meth_states.txt",sep=" ")
rownames(all_CpG_meth) = metadata[which(!is.na(metadata$WGBS)),]$Sample
all_CpG_meth = all_CpG_meth[,2:5]/2
colnames(all_CpG_meth) = meth_states

# TEs
# Number of CpGs overlapping TEs in each methylation state per sample
TE_CpG_meth = read.table("WGBS/CpG_TE_Meth_states.txt",sep=" ")
rownames(TE_CpG_meth) = metadata[which(!is.na(metadata$WGBS)),]$Sample
TE_CpG_meth = TE_CpG_meth[,2:5]/2
colnames(TE_CpG_meth) = meth_states
TE_CpG_meth[is.na(TE_CpG_meth)] = 0

# RefSeq genic features
# Number of CpGs overlapping each RefSeq feature in each methylation state per sample
feature_CpG_meth = read.table("WGBS/feature_CpG_Meth_states.txt",sep="\t",col.names=c("State","Cohort","Sample","CpGs"))
feature_CpG_meth$Sample = mapvalues(feature_CpG_meth$Sample,from = seq(5,41,1), to = as.vector(metadata[which(!is.na(metadata$WGBS)),]$Sample))
feature_CpG_meth$Cohort = gsub("refseq_", "", feature_CpG_meth$Cohort)

# Include all sample x feature x state combinations
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

# Add a column listing the total CpGs per feature and sample
WGBS_state_proportion = ddply(WGBS_state_proportion,.(Sample,Cohort),transform,Total=sum(CpGs))

# DHS

# Genome/TEs
# Number and total width of DHS peaks per sample, overall and overlapping TEs
# Plus number of peaks whose summit overlaps a TE per sample
DNase_stats = read.table("DNase/DNase_stats.txt",sep='\t',col.names=c("Peaks","Sample","Total_width","Peaks_in_TE"))
test = read.table("DNase/rmsk_TEother_merge_DNase_contribution.txt",sep='\t',col.names=c("Sample","Total_width_in_TE"))
DNase_stats = merge(DNase_stats,test,by=c("Sample"))[,c(1:2,4,3,5)]
DNase_stats$Summit_in_TE = read.table("DNase/true_summit/rmsk_TEother_DNase_summit_stats.txt",sep='\t')$V1
rm(test)

# Number of bases of overlap between DHS peaks and each RefSeq genic feature, by sample
DNase_features = read.table("DNase/refseq_features_DNase.txt",sep='\t',col.names=c("Sample","Cohort","Bases"))
DNase_features$Cohort = gsub("refseq_", "", DNase_features$Cohort)

# Add the total width of the feature in that sample
DNase_features$Total = apply(DNase_features,1,function(x) ifelse(x[1] %in% as.vector(metadata[which(metadata$chrY == "Yes"),]$Sample),
                                                                 feature_overlap[which(feature_overlap$Cohort == x[2]),]$Genome,
                                                                 feature_overlap[which(feature_overlap$Cohort == x[2]),]$Genome_noY))

# Combine into one dataframe
DNase_stats_long = melt(DNase_stats[,c("Sample","Total_width","Total_width_in_TE")],
                        id.var=("Sample"),variable.name="Cohort",value.name="Bases")

## Add the total width of the genome/all TEs in that sample
DNase_stats_long$Total = ifelse(DNase_stats_long$Cohort == "Total_width",
                                ifelse(DNase_stats_long$Sample %in% as.vector(metadata[which(metadata$chrY == "Yes"),]$Sample),GENOME_WIDTH,GENOME_WIDTH_noY),
                                ifelse(DNase_stats_long$Sample %in% as.vector(metadata[which(metadata$chrY == "Yes"),]$Sample),MERGED_TE_WIDTH,MERGED_TE_WIDTH_noY))

DNase_stats_long = rbind(DNase_stats_long,DNase_features)
DNase_stats_long$State = rep("DNase",dim(DNase_stats_long)[1])

# H3K27ac 

# Genome/TEs
# Number and total width of H3K27ac peaks per sample, overall and overlapping TEs
# Plus number of peaks whose summit overlaps a TE per sample
H3K27ac_stats = read.table("H3K27ac/H3K27ac_stats.txt",sep='\t',header=TRUE)
test = read.table("H3K27ac/rmsk_TEother_merge_H3K27ac_contribution.txt",sep='\t',col.names=c("Sample","Total_width_in_TE"))
H3K27ac_stats = merge(H3K27ac_stats,test,by=c("Sample"))
H3K27ac_stats$Summit_in_TE = read.table("H3K27ac/true_summit/rmsk_TEother_H3K27ac_summit_stats.txt",sep='\t')$V1
rm(test)

# Number of bases of overlap between H3K27ac peaks and each RefSeq genic feature, by sample
H3K27ac_features = read.table("H3K27ac/refseq_features_H3K27ac.txt",sep='\t',col.names=c("Sample","Cohort","Bases"))
H3K27ac_features$Cohort = gsub("refseq_", "", H3K27ac_features$Cohort)

# Add the total width of the feature in that sample
H3K27ac_features$Total = apply(H3K27ac_features,1,function(x) ifelse(x[1] %in% as.vector(metadata[which(metadata$chrY == "Yes"),]$Sample),
                                                                     feature_overlap[which(feature_overlap$Cohort == x[2]),]$Genome,
                                                                     feature_overlap[which(feature_overlap$Cohort == x[2]),]$Genome_noY))

# Combine into one dataframe
H3K27ac_stats_long = melt(H3K27ac_stats[,c("Sample","Total_width","Total_width_in_TE")],
                          id.var=("Sample"),variable.name="Cohort",value.name="Bases")

## Add the total width of the genome/all TEs in that sample
H3K27ac_stats_long$Total = ifelse(H3K27ac_stats_long$Cohort == "Total_width",
                                  ifelse(H3K27ac_stats_long$Sample %in% as.vector(metadata[which(metadata$chrY == "Yes"),]$Sample),GENOME_WIDTH,GENOME_WIDTH_noY),
                                  ifelse(H3K27ac_stats_long$Sample %in% as.vector(metadata[which(metadata$chrY == "Yes"),]$Sample),MERGED_TE_WIDTH,MERGED_TE_WIDTH_noY))

H3K27ac_stats_long = rbind(H3K27ac_stats_long,H3K27ac_features)
H3K27ac_stats_long$State = rep("H3K27ac",dim(H3K27ac_stats_long)[1])

# Combine matrices of state x sample x feature for all techniques
test = WGBS_state_proportion
colnames(test)[3] = "Bases"
all_state_proportion = rbind(chromHMM_state_proportion,test,DNase_stats_long,H3K27ac_stats_long)
rm(test)

all_state_proportion$Cohort = gsub("Total_width_in_TE","TE",all_state_proportion$Cohort)
all_state_proportion$Cohort = gsub("Total_width","Genome",all_state_proportion$Cohort)
all_state_proportion = all_state_proportion[which(!(all_state_proportion$Cohort %in% c("genes","genes_nc","genes_pc"))),]
all_state_proportion$Mark = factor(c(rep("chromHMM",38100),rep("WGBS",2960),rep("DNase",1060),rep("H3K27ac",1960)),levels=c("chromHMM","WGBS","DNase","H3K27ac"))

# Proportion of feature annotated with each state per sample
all_state_proportion$Proportion = all_state_proportion$Bases/all_state_proportion$Total

# Add metadata
all_state_proportion = merge(all_state_proportion,metadata[,c(1,4:9)],by="Sample")

# Distinguish coding/non-coding features
all_state_proportion = split_coding(all_state_proportion,c("TE","Genome"))

rm(list=c("mnemonics_states_features","feature_CpG_meth","DNase_features","H3K27ac_features","chromHMM_state_proportion","WGBS_state_proportion","DNase_stats_long","H3K27ac_stats_long"))


# Proportion of feature annotated with each state, across samples
## Denominator: chromHMM: total bases annotated; WGBS: total CpGs; DNase/H3K27ac: total feature width
combined_proportion = merge(ddply(all_state_proportion,.(Mark,State,Cohort,Feature,Coding),summarise,Bases=sum(Bases)),
                               ddply(unique(all_state_proportion[,c("Mark","Sample","Cohort","Total")]),.(Mark,Cohort),summarise,Total=sum(Total)),by=c("Mark","Cohort"),all.x=TRUE)
combined_proportion$Proportion = combined_proportion$Bases/combined_proportion$Total
combined_proportion$State = factor(combined_proportion$State,levels=names(all_state_labels)[1:21])


# Proportion of each state that overlap TEs, across samples
## Calculate total bases/CpGs in state for entire genome and overlapping TEs, across samples
contribution = ddply(dcast(all_state_proportion[which(all_state_proportion$Cohort %in% c("Genome","TE")),1:4],Sample+State~Cohort,value.var="Bases"),.(State),summarise,Genome=sum(Genome),TE=sum(TE))
contribution$State = factor(contribution$State,levels=c(levels(contribution$State),"Bases","CpGs"))

## Total bases within TEs (annotated with chromHMM)
contribution[22,] = list("Bases",sum(contribution[1:15,]$Genome),sum(contribution[1:15,]$TE))

## Total CpGs within TEs
contribution[23,] = list("CpGs",ALL_CPGS,TE_CPGS)

# Calculate proportion
contribution$Proportion = contribution$TE/contribution$Genome

# Proportion of each composite chromHMM state that overlap TEs, across samples
# Includes Roadmap active/inactive states as well as custom chromHMM categories
contribution_composite = as.data.frame(rbind(colSums(contribution[1:8,2:3]),colSums(contribution[9:15,2:3]),
                                             colSums(contribution[c(1:3,6:7),2:3]),colSums(contribution[4:5,2:3]),
                                             colSums(contribution[c(9,13:14),2:3]),colSums(contribution[10:12,2:3])))
rownames(contribution_composite) = c("Active","Inactive","Active Regulatory","Transcribed","Repressed","Poised Regulatory")
contribution_composite$Proportion = contribution_composite$TE/contribution_composite$Genome


# Proportion of each state that overlaps TEs, by sample
by_sample_all = dcast(all_state_proportion[which(all_state_proportion$Cohort %in% c("Genome","TE")),1:4],Sample+State~Cohort,value.var="Bases")
by_sample_all$Proportion = by_sample_all$TE/by_sample_all$Genome

## Add metadata
by_sample_all = merge(by_sample_all,metadata[,c(1,4:9)],by="Sample")

## Add total proportion of state within TEs across samples
by_sample_all = merge(by_sample_all,contribution[,c("State","Proportion")],all.x=TRUE,by="State")
colnames(by_sample_all)[c(5,12)] = c("Proportion","Contribution")