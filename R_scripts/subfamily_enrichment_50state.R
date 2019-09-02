# Calculates the log odds ratio enrichment of each TE subfamily in each 50-state model chromHMM state in each sample,
# for 7 samples with individual 50-state chromHMM models.
# Includes thresholds for considering a subfamily x state x sample combination

## subfamily_50state_sample - LOR enrichment for each subfamily x state x sample combination,
## the number of subfamily members overall and in each state per sample, and the corresponding 18-state model state

## subfamily_50state_sample_filter - LOR enrichment for each subfamily x state x sample combination, 
## limited to those that pass thresholds for the number of members overall and in the state

# Length of each state in each subfamily in each sample (ijk)
subfamily_50state_sample = lapply(list.files(path="chromHMM/subfamily/50state/",pattern=".txt",full.names = TRUE),
                                  function(x) read.table(x,sep='\t',col.names=c("subfamily","State","Sample","Length_ijk")))
subfamily_50state_sample = ldply(subfamily_50state_sample)

## Add missing subfamily x state x sample combinations
subfamily_50state_sample_expand = expand.grid(subfamily = rmsk_TE_subfamily$subfamily,
                                              Sample = unique(subfamily_50state_sample$Sample),
                                              State = unique(subfamily_50state_sample$State))
subfamily_50state_sample = merge(subfamily_50state_sample,subfamily_50state_sample_expand,by=c("subfamily","Sample","State"),all.y=TRUE)
subfamily_50state_sample[is.na(subfamily_50state_sample)] = 0
rm(subfamily_50state_sample_expand)

# Length of each state per sample (jk)
sample_50state_files = list.files(path="chromHMM/genome/state50/",pattern="_50_segments_genome.txt")
sample_50state = lapply(sample_50state_files,function(x) read.table(paste("chromHMM/genome/state50/",x,sep=""),sep='\t',col.names=c("State","Length_jk")))
names(sample_50state) = gsub("_50_segments_genome.txt","",sample_50state_files)
sample_50state = ldply(sample_50state,.id="Sample")
rm(sample_50state_files)

# Combine Length ijk and Length jk
subfamily_50state_sample = merge(subfamily_50state_sample,sample_50state,by=c("Sample","State"),all.x=TRUE)
rm(sample_50state)

# Length of subfamily per sample (ik), considering whether the sample includes chrY
subfamily_50state_sample$Length_ik = ifelse(metadata[match(subfamily_50state_sample$Sample,metadata$Sample),]$chrY == "Yes",
                                          rmsk_TE_subfamily[match(subfamily_50state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Total_length,
                                          rmsk_TE_subfamily[match(subfamily_50state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Total_length_noY)

# Length of genome per sample (k), considering whether the sample includes chrY
subfamily_50state_sample$Length_k = ifelse(metadata[match(subfamily_50state_sample$Sample,metadata$Sample),]$chrY == "Yes",GENOME_WIDTH,GENOME_WIDTH_noY)

# Log odds ratio (LOR) enrichment for each subfamily x state x sample combination
# With pseudocount of 1e-20 to avoid undefined values
subfamily_50state_sample$Enrichment = log2(1e-20+((subfamily_50state_sample$Length_ijk/subfamily_50state_sample$Length_jk)/
                                                    (subfamily_50state_sample$Length_ik/subfamily_50state_sample$Length_k)))

# Add class and family information for each subfamily
subfamily_50state_sample = merge(subfamily_50state_sample,rmsk_TE_subfamily[,c("subfamily","family","class_update")],by="subfamily",all.x=TRUE)

# Number of subfamily members in each chromHMM state per sample
## Number of subfamily members where 1) the TE overlaps the center of a 200bp chromHMM bin annotated with the state
## or 2) the state covers the majority of the TE, for TEs that do not overlap the center of any 200bp chromHMM bin,
## by chromHMM state and sample
subfamily_50state_members = lapply(list.files(path="/scratch/ecp/TE_landscape/state50/",pattern="subfamily_state_sample_",full.names = TRUE),
                                   function(x)  read.table(x,sep='\t',col.names = c("subfamily","State","Sample","Members")))
subfamily_50state_members = ldply(subfamily_50state_members)

## Combine the LOR enrichment and number of members in the state for each subfamily x state x sample combination
subfamily_50state_sample = merge(subfamily_50state_sample,subfamily_50state_members,by=c("subfamily","State","Sample"),all.x=TRUE)
subfamily_50state_sample[which(is.na(subfamily_50state_sample$Members)),]$Members = 0
rm(subfamily_50state_members)

# Number of subfamily members per sample, considering whether the sample includes chrY
subfamily_50state_sample$Count = ifelse(metadata[match(subfamily_50state_sample$Sample,metadata$Sample),]$chrY == "Yes",
                                      rmsk_TE_subfamily[match(subfamily_50state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Count,
                                      rmsk_TE_subfamily[match(subfamily_50state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Count_noY)

# For each 50-state model, add the corresponding 18-state model state (from REMC Supplementary Fig. 4)
subfamily_50state_sample = merge(subfamily_50state_sample,state50_state18,by.x=c("Sample","State"),by.y=c("Sample","State50"))
subfamily_50state_sample$State = factor(subfamily_50state_sample$State,levels=unique(state50_state18$State50))
subfamily_50state_sample$State18 = factor(subfamily_50state_sample$State18,levels=names(chromHMM_states_18))

# Filter LOR enrichments to only those where the subfamily has >10 members in the state and >30 overall
subfamily_50state_sample_filter = subfamily_50state_sample[which(subfamily_50state_sample$Members > THRESHOLD_IJK_MEMBER &
                                                                   subfamily_50state_sample$Count > THRESHOLD_IK_MEMBER),]