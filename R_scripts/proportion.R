# Proportion of feature in each state
# See 4/18/2016, 4/20/2016, 4/25/2016, 4/26/2016, 4/27/2016, 5/10/2016, 5/25/2016, 6/27/2016, 9/8/2016, 9/25/2016, 9/27/2016, 2/3/2017, 2/10/2017, 3/8/2017, 5/8/2017, 5/24/2017, 5/30/2017, 6/5/2017, 7/4/2017, 7/24/2017, 8/1/2017, 8/3/2017, 8/7/2017

# chromHMM proportion 

# Genome
# Number of bases in each state in each sample overall
mnemonics_states_genome = read.table("chromHMM/genome/mnemonics_state.txt",sep='\t',header=TRUE,row.names=1)
mnemonics_states_genome$Total = rowSums(mnemonics_states_genome)

# Normalized by total annotated
mnemonics_states_genome_normalized = adply(mnemonics_states_genome[2:16,],1,function(x) x/mnemonics_states_genome[1,])
rownames(mnemonics_states_genome_normalized) = chromHMM_states
mnemonics_states_genome_normalized$Mean = rowMeans(mnemonics_states_genome_normalized[,1:127])
mnemonics_states_genome_normalized$SD = apply(mnemonics_states_genome_normalized[,1:127],1,sd)

# TEs
# Number of bases in each state in each sample in merged TEs
mnemonics_states_TE = read.table("chromHMM/mnemonics_TEother_merge_states.txt",sep='\t',header=TRUE,row.names=1)
mnemonics_states_TE[is.na(mnemonics_states_TE)] = 0
mnemonics_states_TE$Total = rowSums(mnemonics_states_TE)

# Normalized by total annotated
mnemonics_states_TE_normalized = adply(mnemonics_states_TE[2:16,],1,function(x) x/mnemonics_states_TE[1,])
rownames(mnemonics_states_TE_normalized) = chromHMM_states
mnemonics_states_TE_normalized$Mean = rowMeans(mnemonics_states_TE_normalized[,1:127])
mnemonics_states_TE_normalized$SD = apply(mnemonics_states_TE_normalized[,1:127],1,sd)

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

# Normalized across all samples
states_features_noTE = merge(aggregate(data=mnemonics_states_features_noTE,Bases~Cohort+State,sum),aggregate(data=mnemonics_states_features_noTE,Bases~Cohort,sum),by=c("Cohort"),all.x=TRUE)
colnames(states_features_noTE)[3:4] = c("Bases_state","Bases_all")
states_features_noTE$Proportion = states_features_noTE$Bases_state/states_features_noTE$Bases_all

# Combining into one dataframe
mnemonics_states_all = as.data.frame(c(mnemonics_states_genome_normalized$Total,mnemonics_states_TE_normalized$Total))
colnames(mnemonics_states_all)[1] = "Proportion"
mnemonics_states_all$State = rep(chromHMM_states,2)
mnemonics_states_all$Cohort = c(rep("Genome",15),rep("TE",15))
mnemonics_states_all = rbind(mnemonics_states_all,states_features_noTE[,c(1:2,5)])

# WGBS proportion 

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

# Proportion of Refseq feature and TE CpGs hypomethylated per sample
CpG_Meth = as.data.frame(rbind(colSums(all_CpG_meth),colSums(TE_CpG_meth)))
CpG_Meth$Cohort = c("Genome","TE")
CpG_Meth = as.data.frame(rbind(CpG_Meth,aggregate(data=feature_CpG_meth[,1:5],.~Cohort,sum)))
rownames(CpG_Meth) = CpG_Meth$Cohort
CpG_Meth = CpG_Meth[,1:4]
CpG_Meth = melt(as.matrix(CpG_Meth/rowSums(CpG_Meth)))
colnames(CpG_Meth) = c("Cohort","State","Proportion")

# DNase proportion

# Number and width of DNase peaks per sample, overall and in TEs 
DNase_stats = read.table("DNase/DNase_stats.txt",sep='\t')
colnames(DNase_stats) = c("Peaks","Sample","Total_width","Peaks_in_TE")
test = read.table("DNase/rmsk_TEother_merge_DNase_contribution.txt",sep='\t')
colnames(test) = c("Sample","Total_width_in_TE")
DNase_stats = merge(DNase_stats,test,by=c("Sample"))[,c(1:2,4,3,5)]
rm(test)

# Number of bases of DNase peak in Refseq genic features, by sample
DNase_features = read.table("DNase/Refseq_features/refseq_features_DNase.txt",sep='\t')
colnames(DNase_features) = c("Sample","Bases","Cohort")
DNase_features$Cohort = gsub("refseq_", "", DNase_features$Cohort)
DNase_features$Cohort = gsub("_merge_noTE_DNase", "", DNase_features$Cohort)
DNase_features$Cohort = gsub("_noTE_DNase", "", DNase_features$Cohort)
DNase_features$Cohort = gsub("coding_exon", "CDS", DNase_features$Cohort)
DNase_features = DNase_features[which(!(DNase_features$Cohort %in% c("genes","genes_nc","genes_pc"))),]

# Proportion of RefSeq features overlapping DNase peaks, averaged
DNase_features$Genome = apply(DNase_features,1,function(x) ifelse(x[1] %in% as.vector(metadata[which(metadata$chrY == "Yes"),]$Sample),
                                                                  feature_overlap[which(feature_overlap$Cohort == x[3]),]$Genome,
                                                                  feature_overlap[which(feature_overlap$Cohort == x[3]),]$Genome_noY))
colnames(DNase_features)[c(2,4)] = c("Total_width_in_feature","Feature_width")
DNase_features$Proportion = DNase_features$Total_width_in_feature/DNase_features$Feature_width

# Proportion of TEs overlapping DNase peaks, across all samples
DNase_proportion = aggregate(data=DNase_features[,2:4],.~Cohort,sum)
DNase_proportion[20,] = c("Genome",sum(as.numeric(DNase_stats$Total_width)),DNASE_TOTAL_WIDTH)
DNase_proportion[21,] = c("TE",sum(as.numeric(DNase_stats$Total_width_in_TE)),DNASE_TE_WIDTH)
DNase_proportion$Proportion = as.numeric(DNase_proportion$Total_width_in_feature)/as.numeric(DNase_proportion$Feature_width)
DNase_proportion$State = rep("DNase",dim(DNase_proportion)[1])

# H3K27ac proportion

# Number and width of H3K27ac peaks per sample, overall and in TEs 
H3K27ac_stats = read.table("H3K27ac/H3K27ac_stats.txt",sep='\t',header=TRUE)
test = read.table("H3K27ac/rmsk_TEother_merge_H3K27ac_contribution.txt",sep='\t')
colnames(test) = c("Sample","Total_width_in_TE")
H3K27ac_stats = merge(H3K27ac_stats,test,by=c("Sample"))
rm(test)

# Number of bases of H3K27ac peak in Refseq genic features, by sample
H3K27ac_features = read.table("H3K27ac/Refseq_features/refseq_features_H3K27ac.txt",sep='\t')
colnames(H3K27ac_features) = c("Sample","Bases","Cohort")
H3K27ac_features$Cohort = gsub("refseq_", "", H3K27ac_features$Cohort)
H3K27ac_features$Cohort = gsub("_merge_noTE_H3K27ac", "", H3K27ac_features$Cohort)
H3K27ac_features$Cohort = gsub("_noTE_H3K27ac", "", H3K27ac_features$Cohort)
H3K27ac_features$Cohort = gsub("coding_exon", "CDS", H3K27ac_features$Cohort)

# Proportion of RefSeq features overlapping H3K27ac peaks, averaged
H3K27ac_features$Genome = apply(H3K27ac_features,1,function(x) ifelse(x[1] %in% as.vector(metadata[which(metadata$chrY == "Yes"),]$Sample),
                                                                  feature_overlap[which(feature_overlap$Cohort == x[3]),]$Genome,
                                                                  feature_overlap[which(feature_overlap$Cohort == x[3]),]$Genome_noY))
colnames(H3K27ac_features)[c(2,4)] = c("Total_width_in_feature","Feature_width")
H3K27ac_features$Proportion = H3K27ac_features$Total_width_in_feature/H3K27ac_features$Feature_width

# Proportion of TEs overlapping H3K27ac peaks, across all samples
H3K27ac_proportion = aggregate(data=H3K27ac_features[,2:4],.~Cohort,sum)
H3K27ac_proportion[20,] = c("Genome",sum(as.numeric(H3K27ac_stats$Total_width)),H3K27AC_TOTAL_WIDTH)
H3K27ac_proportion[21,] = c("TE",sum(as.numeric(H3K27ac_stats$Total_width_in_TE)),H3K27AC_TE_WIDTH)
H3K27ac_proportion$Proportion = as.numeric(H3K27ac_proportion$Total_width_in_feature)/as.numeric(H3K27ac_proportion$Feature_width)
H3K27ac_proportion$State = rep("H3K27ac",dim(H3K27ac_proportion)[1])

# Combine all
combined_proportion = rbind(mnemonics_states_all,CpG_Meth,DNase_proportion[,c(1,4:5)],H3K27ac_proportion[,c(1,4:5)])
combined_proportion$State = factor(combined_proportion$State,levels=names(all_state_labels)[1:21])
combined_proportion$Group = factor(c(rep("chromHMM",315),rep("WGBS",80),rep("DNase",21),rep("H3K27ac",21)),levels=c("chromHMM","WGBS","DNase","H3K27ac"))
combined_proportion = split_coding(combined_proportion,c("TE","Genome","genome_noTE"))