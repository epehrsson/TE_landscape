# Calculates the log odds ratio enrichment of each TE subfamily in each epigenetic state in each sample
# and the proportion of each state within each subfamily by sample, then creates tables
# of the number of samples each subfamily is enriched LOR > 1.5 or represents >1% of an epigenetic state.
# Includes thresholds for considering a subfamily x state x sample combination

## TE_meth_subfamily - number of members of each subfamily in each methylation state per sample
## subfamily_state_sample_combined - LOR enrichment for each subfamily x state x sample combination, 
## the proportion of each state in each subfamily by sample, and the number of subfamily members overall and in each state per sample

## subfamily_state_sample_filter - LOR enrichment for each subfamily x state x sample combination, 
## limited to those that pass thresholds for the number of members overall and in the state

## subfamily_state_sample_counts - number/proportion of samples with LOR enrichment > 1.5 per subfamily x state, 
## including only those that pass subfamily member thresholds

## subfamily_state_sample_counts_combine - number of subfamilies with at least one LOR enrichment > 1.5 per state, overall and by class
## including only those that pass subfamily member thresholds

## subfamily_state_sample_counts_pc - number/proportion of samples in which each subfamily represents >1% of each state
## subfamily_state_sample_counts_pc_combine - number of subfamilies that represent >1% of each state at least once, overall and by class


# Length of each state in each subfamily in each sample (ijk)

## chromHMM: length of overlap (bp) between subfamily and chromHMM state per sample
subfamily_chromHMM_sample = read.table("chromHMM/subfamily/subfamily_state_sample.txt",sep='\t',col.names=c("subfamily","State","Sample","Length_ijk"))
### Add missing subfamily x state x sample combinations
subfamily_chromHMM_sample_expand = expand.grid(subfamily = rmsk_TE_subfamily$subfamily,Sample = metadata$Sample,State = chromHMM_states)
subfamily_chromHMM_sample = merge(subfamily_chromHMM_sample,subfamily_chromHMM_sample_expand,by=c("subfamily","Sample","State"),all.y=TRUE)
subfamily_chromHMM_sample[is.na(subfamily_chromHMM_sample)] = 0
rm(subfamily_chromHMM_sample_expand)

## DHS: number of unique peaks whose summit overlaps each subfamily per sample
subfamily_DNase_sample = read.table("DNase/true_summit/subfamily_DNase_sample_summit.txt",sep='\t',col.names=c("Sample","subfamily","Length_ijk"))
subfamily_DNase_sample$State = rep("DNase",dim(subfamily_DNase_sample)[1])

## H3K27ac: number of unique peaks whose summit overlaps each subfamily per sample
subfamily_H3K27ac_sample = read.table("H3K27ac/true_summit/subfamily_H3K27ac_sample_summit.txt",sep='\t',col.names=c("Sample","subfamily","Length_ijk"))
subfamily_H3K27ac_sample$State = rep("H3K27ac",dim(subfamily_H3K27ac_sample)[1])

## Combine into a single dataframe
subfamily_state_sample = rbind(subfamily_chromHMM_sample,subfamily_DNase_sample,subfamily_H3K27ac_sample)
rm(list=c("subfamily_chromHMM_sample","subfamily_DNase_sample","subfamily_H3K27ac_sample"))

# Length of each state per sample (jk)

## chromHMM: length of each chromHMM state (bp) per sample
chromHMM_state_sample = mnemonics_states_genome[,1:3]
colnames(chromHMM_state_sample)[3] = "Length_jk"

## DHS: number of unique peaks per sample
DNase_peaks_sample = DNase_stats[,c("Sample","Peaks")]
DNase_peaks_sample$State = rep("DNase",dim(DNase_peaks_sample)[1])
colnames(DNase_peaks_sample)[2] = "Length_jk"

## H3K27ac: number of unique peaks per sample
H3K27ac_peaks_sample = H3K27ac_stats[,c("Sample","Peaks")]
H3K27ac_peaks_sample$State = rep("H3K27ac",dim(H3K27ac_peaks_sample)[1])
colnames(H3K27ac_peaks_sample)[2] = "Length_jk"

## Combine into a single dataframe
peaks_sample = rbind(chromHMM_state_sample,DNase_peaks_sample,H3K27ac_peaks_sample)
rm(list=c("chromHMM_state_sample","DNase_peaks_sample","H3K27ac_peaks_sample"))

# Combine Length ijk and Length jk
subfamily_state_sample = merge(subfamily_state_sample,peaks_sample,by=c("Sample","State"),all.x=TRUE)
rm(peaks_sample)

# Length of subfamily per sample (ik), considering whether the sample includes chrY
subfamily_state_sample$Length_ik = ifelse(metadata[match(subfamily_state_sample$Sample,metadata$Sample),]$chrY == "Yes",
                                          rmsk_TE_subfamily[match(subfamily_state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Total_length,
                                          rmsk_TE_subfamily[match(subfamily_state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Total_length_noY)

# Length of genome per sample (k), considering whether the sample includes chrY
subfamily_state_sample$Length_k = ifelse(metadata[match(subfamily_state_sample$Sample,metadata$Sample),]$chrY == "Yes",GENOME_WIDTH,GENOME_WIDTH_noY)

# Log odds ratio (LOR) enrichment for each subfamily x state x sample combination
# With pseudocount of 1e-20 to avoid undefined values
subfamily_state_sample$Enrichment = log2(1e-20+((subfamily_state_sample$Length_ijk/subfamily_state_sample$Length_jk)/(subfamily_state_sample$Length_ik/subfamily_state_sample$Length_k)))

# Add class and family information for each subfamily
subfamily_state_sample = merge(subfamily_state_sample,rmsk_TE_subfamily[,c("subfamily","family","class_update")],by="subfamily",all.x=TRUE)

# Add metadata for each sample
subfamily_state_sample = merge(subfamily_state_sample,metadata[,c(1,4:9)],by=c("Sample"),all.x=TRUE)

# Proportion of each epigenetic state in each subfamily per sample
subfamily_state_sample$Length_percent_jk = subfamily_state_sample$Length_ijk/subfamily_state_sample$Length_jk

# Number of subfamily members in each epigenetic state per sample

## chromHMM: number of subfamily members where 1) the TE overlaps the center of a 200bp chromHMM bin annotated with the state
## or 2) the state covers the majority of the TE, for TEs that do not overlap the center of any 200bp chromHMM bin,
## by chromHMM state and sample
subfamily_chromHMM_sample_members = read.table("chromHMM/subfamily/subfamily_state_sample_summit.txt",sep='\t',col.names=c("subfamily","Sample","State","Members"))

## DHS: number of subfamily members overlapping the summit of a DHS peak per sample (column)
TE_DNase_peaks_members = melt(aggregate(data=TE_DNase_peaks[,c(4,8:60)],.~subfamily,function(x) sum(as.numeric(x) > 0)),id.vars=c("subfamily"),
                              variable.name="Sample",value.name="Members")
TE_DNase_peaks_members$State = rep("DNase",dim(TE_DNase_peaks_members)[1])

## H3K27ac: number of subfamily members overlapping the summit of an H3K27ac peak per sample (column)
TE_H3K27ac_peaks_members = melt(aggregate(data=TE_H3K27ac_peaks[,c(4,8:105)],.~subfamily,function(x) sum(as.numeric(x) > 0)),id.vars=c("subfamily"),
                                variable.name="Sample",value.name="Members")
TE_H3K27ac_peaks_members$State = rep("H3K27ac",dim(TE_H3K27ac_peaks_members)[1])

## Combine into a single dataframe
subfamily_state_sample_members = rbind(subfamily_chromHMM_sample_members,TE_DNase_peaks_members,TE_H3K27ac_peaks_members)
rm(list=c("subfamily_chromHMM_sample_members","TE_DNase_peaks_members","TE_H3K27ac_peaks_members"))

## Combine the LOR enrichment and number of members in the state for each subfamily x state x sample combination
subfamily_state_sample = merge(subfamily_state_sample,subfamily_state_sample_members,by=c("subfamily","State","Sample"),all.x=TRUE)
subfamily_state_sample[which(is.na(subfamily_state_sample$Members)),]$Members = 0
rm(subfamily_state_sample_members)

# Proportion of all subfamily members in each epigenetic state per sample, considering whether the sample includes chrY
subfamily_state_sample$Count = ifelse(metadata[match(subfamily_state_sample$Sample,metadata$Sample),]$chrY == "Yes",
                                                                                       rmsk_TE_subfamily[match(subfamily_state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Count,
                                                                                       rmsk_TE_subfamily[match(subfamily_state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Count_noY)
subfamily_state_sample$Percent = subfamily_state_sample$Members/subfamily_state_sample$Count


# Number of unique CpGs overlapping each subfamily in each methylation state per sample (ijk)
subfamily_CpG_meth = read.table("WGBS/subfamily_CpG_Meth_states.txt",sep='\t',col.names=c("Sample",meth_states,"subfamily"))
## Add sample names
subfamily_CpG_meth$Sample = mapvalues(subfamily_CpG_meth$Sample,seq(4,40,1),as.vector(metadata[which(!is.na(metadata$WGBS)),]$Sample))
subfamily_CpG_meth[is.na(subfamily_CpG_meth)] = 0
subfamily_CpG_meth[,meth_states] = subfamily_CpG_meth[,meth_states]/2
subfamily_CpG_meth = melt(subfamily_CpG_meth,id.vars=c("Sample","subfamily"),
                          variable.name="State",value.name="CpG_ijk")

# Number of CpGs per subfamily (ik)
subfamily_CpG_meth = merge(subfamily_CpG_meth,rmsk_TE_subfamily[,c("subfamily","family","class_update","CpGs")],by=c("subfamily"),all.x=TRUE)
colnames(subfamily_CpG_meth)[7] = "CpG_ik"

# Number of CpGs in each methylation state per sample (jk)
subfamily_CpG_meth$CpG_jk = apply(subfamily_CpG_meth,1,function(x) all_CpG_meth[x[2],x[3]])

# Number of CpGs in genome (k)
subfamily_CpG_meth$CpG_k = rep(ALL_CPGS,dim(subfamily_CpG_meth)[1])

# Log odds ratio (LOR) enrichment for each subfamily x methylation state x sample combination
# With pseudocount of 1e-20 to avoid undefined values
subfamily_CpG_meth$Enrichment = log2(1e-20+((subfamily_CpG_meth$CpG_ijk/subfamily_CpG_meth$CpG_ik)/(subfamily_CpG_meth$CpG_jk/subfamily_CpG_meth$CpG_k)))

# Proportion of CpGs in each methylation state in each subfamily per sample
subfamily_CpG_meth$CpG_ijk_jk = subfamily_CpG_meth$CpG_ijk/subfamily_CpG_meth$CpG_jk

# Add metadata for each sample
subfamily_CpG_meth = merge(subfamily_CpG_meth,metadata[,c(1,4:9)],by=c("Sample"),all.x=TRUE)

# Number of subfamily members in each methylation state per sample, based on TE average methylation
## Average methylation of each TE in each sample
TE_meth_average_long = melt(TE_meth_average[,c(4,6,8:44,54)],id.vars=c("subfamily","family","class_update"),
                            variable.name="Sample",value.name="Methylation")
## Assign each TE a methylation state based on average methylation, 
## then count the number of TEs per subfamily and sample in each methylation state
TE_meth_subfamily = ddply(TE_meth_average_long,.(subfamily,family,class_update,Sample),summarise,Hypomethylated=sum(na.omit(Methylation) < 0.3),
                          Intermediate=sum(na.omit(Methylation) >= 0.3 & na.omit(Methylation) <= 0.7),Hypermethylated=sum(na.omit(Methylation) > 0.7),
                          Missing=sum(is.na(Methylation)))
TE_meth_subfamily = melt(TE_meth_subfamily,id.vars=c("subfamily","family","class_update","Sample"),
                         variable.name="State",value.name="Members")
rm(TE_meth_average_long)

## Combine the LOR enrichment and number of members in the state for each subfamily x methylation state x sample combination
subfamily_CpG_meth = merge(subfamily_CpG_meth,TE_meth_subfamily[,c("subfamily","State","Sample","Members")],by=c("subfamily","Sample","State"))

# Proportion of all subfamily members in each methylation state per sample
# Including only subfamily members that overlap at least one CpG
subfamily_CpG_meth$Count = rmsk_TE_subfamily[match(subfamily_CpG_meth$subfamily,rmsk_TE_subfamily$subfamily),]$Count_CpGs
subfamily_CpG_meth$Percent = subfamily_CpG_meth$Members/subfamily_CpG_meth$Count


# Combine LOR enrichment matrices for all four epigenetic marks
columns = c("subfamily","family","class_update","State","Sample",sample_categories,enrichment_names)
subfamily_state_sample_combined = rbind(subfamily_state_sample[,columns],
                                        rename(subfamily_CpG_meth,c("CpG_ijk"="Length_ijk","CpG_ik"="Length_ik","CpG_jk"="Length_jk","CpG_k"="Length_k",
                                                                    "CpG_ijk_jk"="Length_percent_jk"))[,columns])

# Filter LOR enrichments to only those where the subfamily has >10 members in the state and >30 overall
# For the methylation states, considers only TEs that overlap at least one CpG
subfamily_state_sample_filter = subfamily_state_sample_combined[which(subfamily_state_sample_combined$Members > THRESHOLD_IJK_MEMBER & subfamily_state_sample_combined$Count > THRESHOLD_IK_MEMBER),]
subfamily_state_sample_filter$State = factor(subfamily_state_sample_filter$State,levels=states[1:21])


# Number of LOR enrichments > 1.5 per subfamily x state, including only those that pass subfamily member thresholds
subfamily_state_sample_counts = ddply(subfamily_state_sample_filter,.(class_update,family,subfamily,State),function(x) sum(x$Enrichment > THRESHOLD_LOR))
## Add missing subfamily x state combinations
subfamily_state_expand = expand.grid(subfamily = levels(subfamily_state_sample$subfamily),State = levels(subfamily_state_sample_combined$State))
subfamily_state_expand = join(subfamily_state_expand,rmsk_TE_subfamily[,c("subfamily","family","class_update")],by=c("subfamily"))
subfamily_state_sample_counts = join(subfamily_state_expand,subfamily_state_sample_counts,by=c("subfamily","family","class_update","State"),type="left")[,c(1,3,4,2,5)]
rm(subfamily_state_expand)
subfamily_state_sample_counts[is.na(subfamily_state_sample_counts)] = 0
subfamily_state_sample_counts$State = factor(subfamily_state_sample_counts$State,levels=states[1:21])

# Proportion of samples in which each TE subfamily is enriched in each epigenetic state
# Considering only samples with data for each epigenetic mark
subfamily_state_sample_counts$Sample.Proportion = ifelse(subfamily_state_sample_counts$State %in% chromHMM_states,subfamily_state_sample_counts$V1/sample_counts["All","chromHMM"],
                                                         ifelse(subfamily_state_sample_counts$State %in% meth_states,subfamily_state_sample_counts$V1/sample_counts["All","WGBS"],
                                                                ifelse(subfamily_state_sample_counts$State == "DNase",subfamily_state_sample_counts$V1/sample_counts["All","DNase"],
                                                                       ifelse(subfamily_state_sample_counts$State == "H3K27ac",subfamily_state_sample_counts$V1/sample_counts["All","H3K27ac"],"NA"))))
subfamily_state_sample_counts$Sample.Proportion = as.numeric(subfamily_state_sample_counts$Sample.Proportion)

# Number of subfamilies with at least one enrichment > 1.5 per state, overall and by class
subfamily_state_sample_counts_combine = merge(aggregate(data=subfamily_state_sample_counts,V1~State,function(x) sum(x > 0)),
                                              dcast(aggregate(data=subfamily_state_sample_counts,V1~State+class_update,function(x) sum(x > 0)),State~class_update,value.var = "V1"),
                                              by=c("State"))
subfamily_state_sample_counts_combine = subfamily_state_sample_counts_combine[match(levels(subfamily_state_sample_counts$State),subfamily_state_sample_counts_combine$State),]


# Number of samples in which each subfamily represents >1% of each state
subfamily_state_sample_counts_pc = ddply(subfamily_state_sample_combined,.(class_update,family,subfamily,State),function(x) sum(x$Length_percent_jk > THRESHOLD_PC))

# Proportion of samples in which each TE subfamily represents >1% of each epigenetic state
# Considering only samples with data for each epigenetic mark
subfamily_state_sample_counts_pc$Sample.Proportion = ifelse(subfamily_state_sample_counts_pc$State %in% chromHMM_states,subfamily_state_sample_counts_pc$V1/sample_counts["All","chromHMM"],
                                                            ifelse(subfamily_state_sample_counts_pc$State %in% meth_states,subfamily_state_sample_counts_pc$V1/sample_counts["All","WGBS"],
                                                                   ifelse(subfamily_state_sample_counts_pc$State == "DNase",subfamily_state_sample_counts_pc$V1/sample_counts["All","DNase"],
                                                                          ifelse(subfamily_state_sample_counts_pc$State == "H3K27ac",subfamily_state_sample_counts_pc$V1/sample_counts["All","H3K27ac"],"NA"))))
subfamily_state_sample_counts_pc$Sample.Proportion = as.numeric(subfamily_state_sample_counts_pc$Sample.Proportion)

# Number of subfamilies that represent >1% of each state at least once, overall and by class
subfamily_state_sample_counts_pc_combine = merge(aggregate(data=subfamily_state_sample_counts_pc,V1~State,function(x) sum(x > 0)),
                                              dcast(aggregate(data=subfamily_state_sample_counts_pc,V1~State+class_update,function(x) sum(x > 0)),State~class_update,value.var = "V1"),
                                              by=c("State"))

rm(list=c("subfamily_state_sample","subfamily_CpG_meth"))
