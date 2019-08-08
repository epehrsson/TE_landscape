# Subfamily enrichment using the 50-state models

# Length of state in subfamily in sample (ijk)
subfamily_50state_sample = lapply(list.files(path="chromHMM/subfamily/50state/",pattern=".txt",full.names = TRUE),
                                  function(x) read.table(x,sep='\t',col.names=c("subfamily","State","Sample","Length_ijk")))
subfamily_50state_sample = ldply(subfamily_50state_sample)

## Add missing combinations
subfamily_50state_sample_expand = expand.grid(subfamily = rmsk_TE_subfamily$subfamily,
                                              Sample = unique(subfamily_50state_sample$Sample),
                                              State = unique(subfamily_50state_sample$State))
subfamily_50state_sample = merge(subfamily_50state_sample,subfamily_50state_sample_expand,by=c("subfamily","Sample","State"),all.y=TRUE)
subfamily_50state_sample[is.na(subfamily_50state_sample)] = 0
rm(subfamily_50state_sample_expand)

# Length of state in sample (jk)
sample_50state_files = list.files(path="chromHMM/genome/state50/",pattern="_50_segments_genome.txt")
sample_50state = lapply(sample_50state_files,function(x) read.table(paste("chromHMM/genome/state50/",x,sep=""),sep='\t',col.names=c("State","Length_jk")))
names(sample_50state) = gsub("_50_segments_genome.txt","",sample_50state_files)
sample_50state = ldply(sample_50state,.id="Sample")
rm(sample_50state_files)

## Join
subfamily_50state_sample = merge(subfamily_50state_sample,sample_50state,by=c("Sample","State"),all.x=TRUE)
rm(sample_50state)

# Length of subfamily (ik)
subfamily_50state_sample$Length_ik = ifelse(metadata[match(subfamily_50state_sample$Sample,metadata$Sample),]$chrY == "Yes",
                                          rmsk_TE_subfamily[match(subfamily_50state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Total_length,
                                          rmsk_TE_subfamily[match(subfamily_50state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Total_length_noY)

# Length of genome (k)
subfamily_50state_sample$Length_k = ifelse(metadata[match(subfamily_50state_sample$Sample,metadata$Sample),]$chrY == "Yes",GENOME_WIDTH,GENOME_WIDTH_noY)

# LOR Enrichment for subfamily x state x sample
subfamily_50state_sample$Enrichment = log2(1e-20+((subfamily_50state_sample$Length_ijk/subfamily_50state_sample$Length_jk)/
                                                    (subfamily_50state_sample$Length_ik/subfamily_50state_sample$Length_k)))

# Add class and family
subfamily_50state_sample = merge(subfamily_50state_sample,rmsk_TE_subfamily[,c("subfamily","family","class_update")],by="subfamily",all.x=TRUE)

# Number of subfamily members in state
subfamily_50state_members = lapply(list.files(path="/scratch/ecp/TE_landscape/state50/",pattern="subfamily_state_sample_",full.names = TRUE),
                                   function(x)  read.table(x,sep='\t',col.names = c("subfamily","State","Sample","Members")))
subfamily_50state_members = ldply(subfamily_50state_members)

## Join
subfamily_50state_sample = merge(subfamily_50state_sample,subfamily_50state_members,by=c("subfamily","State","Sample"),all.x=TRUE)
subfamily_50state_sample[which(is.na(subfamily_50state_sample$Members)),]$Members = 0
rm(subfamily_50state_members)

# Number of subfamily members total
subfamily_50state_sample$Count = ifelse(metadata[match(subfamily_50state_sample$Sample,metadata$Sample),]$chrY == "Yes",
                                      rmsk_TE_subfamily[match(subfamily_50state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Count,
                                      rmsk_TE_subfamily[match(subfamily_50state_sample$subfamily,rmsk_TE_subfamily$subfamily),]$Count_noY)

# Add corresponding 18-state model state
subfamily_50state_sample = merge(subfamily_50state_sample,state50_state18,by.x=c("Sample","State"),by.y=c("Sample","State50"))
subfamily_50state_sample$State = factor(subfamily_50state_sample$State,levels=unique(state50_state18$State50))
subfamily_50state_sample$State18 = factor(subfamily_50state_sample$State18,levels=names(chromHMM_states_18))

# Filter matrix
subfamily_50state_sample_filter = subfamily_50state_sample[which(subfamily_50state_sample$Members > THRESHOLD_IJK_MEMBER &
                                                                   subfamily_50state_sample$Count > THRESHOLD_IK_MEMBER),]