# Generates dataframes with the proportion of each shared subfamily hypomethylated, in each chromHMM state, or in an active regulatory state
# in each species in each human/mouse sample pair, restricted to those with >30 members in both species

## hg19_mm10_subfamily_count - For subfamilies in both hg19 and mm10, the number of members in both species, overall and overlapping a CpG
## hg19_mm10_TE_WGBS_subfamily_hypo_paired - Proportion of members hypomethylated in each shared subfamily in each species in each human/mouse sample pair
## hg19_mm10_chromHMM_subfamily - Proportion of each shared subfamily in each chromHMM state in each species in each human/mouse sample pair
## hg19_mm10_chromHMM_subfamily_active - Proportion of each shared subfamily in an active regulatory chromHMM state in each species in each human/mouse sample pair


# Number of members per mm10 subfamily, overall and overlapping a CpG
mm10_subfamily_count = merge(ddply(mm10_rmsk_TE,.(subfamily),summarise,Count=length(subfamily)),
                             ddply(mm10_rmsk_TE_WGBS,.(subfamily),summarise,Count_CpGs=length(subfamily)),by="subfamily",all.x=TRUE)
mm10_subfamily_count[is.na(mm10_subfamily_count)] = 0

# For subfamilies in both hg19 and mm10, the number of members in both species, overall and overlapping a CpG
hg19_mm10_subfamily_count = merge(rmsk_TE_subfamily[,c("subfamily","family","class_update","Count","Count_CpGs")],mm10_subfamily_count,by="subfamily")


# Proportion of members hypomethylated in each subfamily in each species in each human/mouse sample pair,
# Restricted to subfamilies with > 30 members with CpGs in both species
## Proportion of members hypomethylated in each mm10 subfamily (row) in each sample (column)
mm10_rmsk_TE_WGBS_subfamily_hypo = aggregate(data=mm10_rmsk_TE_WGBS[,c("subfamily",as.vector(na.omit(human_mouse_samples)$Mouse_sample_WGBS))],
                                             .~subfamily,function(x) sum(na.omit(x) < 0.3)/length(x),na.action=na.pass)

## Proportion of members hypomethylated in each hg19 subfamily (row) in each sample (column)
hg19_rmsk_TE_WGBS_subfamily_hypo = merge(TE_meth_subfamily[which(TE_meth_subfamily$State == "Hypomethylated" & 
                                                                   TE_meth_subfamily$Sample %in% as.vector(na.omit(human_mouse_samples)$Human_sample)),
                                                     c("subfamily","family","class_update","Sample","Members")],rmsk_TE_subfamily[,c("subfamily","Count_CpGs")],by="subfamily")
hg19_rmsk_TE_WGBS_subfamily_hypo$Percent = hg19_rmsk_TE_WGBS_subfamily_hypo$Members/hg19_rmsk_TE_WGBS_subfamily_hypo$Count_CpGs
hg19_rmsk_TE_WGBS_subfamily_hypo = dcast(hg19_rmsk_TE_WGBS_subfamily_hypo,subfamily+family+class_update~Sample,value.var="Percent")

## Combine
hg19_mm10_TE_WGBS_subfamily_hypo = merge(hg19_rmsk_TE_WGBS_subfamily_hypo,mm10_rmsk_TE_WGBS_subfamily_hypo,by="subfamily")
rm(list=c("hg19_rmsk_TE_WGBS_subfamily_hypo","mm10_rmsk_TE_WGBS_subfamily_hypo"))

## Restrict to subfamilies with >30 TEs with CpGs in both human and mouse
hg19_mm10_TE_WGBS_subfamily_hypo = hg19_mm10_TE_WGBS_subfamily_hypo[which(hg19_mm10_TE_WGBS_subfamily_hypo$subfamily %in% 
                                                                            as.vector(hg19_mm10_subfamily_count[which(
                                                                              hg19_mm10_subfamily_count$Count_CpGs.x > THRESHOLD_IK_MEMBER & 
                                                                                hg19_mm10_subfamily_count$Count_CpGs.y > THRESHOLD_IK_MEMBER),]$subfamily)),]

## Reformat and restrict to human/mouse sample pairs
hg19_mm10_TE_WGBS_subfamily_hypo_paired = melt(hg19_mm10_TE_WGBS_subfamily_hypo,
                                               id.vars=c("subfamily","family","class_update",as.vector(na.omit(human_mouse_samples)$Human_sample)),
                                               variable.name="Mouse_sample_WGBS", value.name="Mouse_hypo")
hg19_mm10_TE_WGBS_subfamily_hypo_paired = melt(hg19_mm10_TE_WGBS_subfamily_hypo_paired,
                                               id.vars=c("subfamily","family","class_update","Mouse_sample_WGBS","Mouse_hypo"),
                                               variable.name="Human_sample",value.name="Human_hypo")
hg19_mm10_TE_WGBS_subfamily_hypo_paired = merge(hg19_mm10_TE_WGBS_subfamily_hypo_paired,human_mouse_samples,by=c("Mouse_sample_WGBS","Human_sample"))

## Calculate the ratio of hypomethylated members in human versus mouse for each subfamily
hg19_mm10_TE_WGBS_subfamily_hypo_paired$Ratio = hg19_mm10_TE_WGBS_subfamily_hypo_paired$Human_hypo/hg19_mm10_TE_WGBS_subfamily_hypo_paired$Mouse_hypo
rm(hg19_mm10_TE_WGBS_subfamily_hypo)


# Proportion of each shared subfamily in each chromHMM state in each species in each human/mouse sample pair
## Proportion of each hg19 subfamily in each chromHMM state in each sample, excluding those with 30 or fewer members
hg19_chromHMM_subfamily = subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% chromHMM_states & 
                                                                  subfamily_state_sample_combined$Sample %in% as.vector(human_mouse_samples$Human_sample) & 
                                                                  subfamily_state_sample_combined$Count > THRESHOLD_IK_MEMBER),
                                                          c("subfamily","family","class_update","Sample","State","Members","Count","Percent")]
colnames(hg19_chromHMM_subfamily)[4:5] = c("Human_sample","Human_state_chromHMM")

## Proportion of each mm10 subfamily in each chromHMM state in each sample, excluding those with 30 or fewer members
### Load number of members in each state per subfamily x sample
mm10_chromHMM_subfamily = read.table("Mouse/chromHMM/mm10_chromHMM_subfamily.txt",sep='\t',
                                     col.names=c("subfamily","Mouse_state_chromHMM","Mouse_sample_chromHMM","Members"))
### Add missing subfamily x state x sample combinations
mm10_chromHMM_subfamily = merge(mm10_chromHMM_subfamily,
                                expand.grid(subfamily=as.vector(unique(mm10_rmsk_TE$subfamily)),
                                            Mouse_state_chromHMM=mouse_chromHMM_states,
                                            Mouse_sample_chromHMM=as.vector(human_mouse_samples$Mouse_sample_chromHMM)),
                                by=c("subfamily","Mouse_state_chromHMM","Mouse_sample_chromHMM"),all.y=TRUE)
mm10_chromHMM_subfamily[is.na(mm10_chromHMM_subfamily)] = 0
### Add number of subfamily members
mm10_chromHMM_subfamily = merge(mm10_chromHMM_subfamily,mm10_subfamily_count[,c("subfamily","Count")],by="subfamily")
### Calculate proprotion of subfamily members in state per sample
mm10_chromHMM_subfamily$Percent = mm10_chromHMM_subfamily$Members/mm10_chromHMM_subfamily$Count
### Restrict to subfamilies with > 30 members 
mm10_chromHMM_subfamily = mm10_chromHMM_subfamily[which(mm10_chromHMM_subfamily$Count > THRESHOLD_IK_MEMBER),]

## Combine, restricting to human/mouse sample pairs, shared subfamilies only
hg19_mm10_chromHMM_subfamily = apply(human_mouse_samples,1,function(x) merge(hg19_chromHMM_subfamily[which(hg19_chromHMM_subfamily$Human_sample == x[1]),],
                                                                             mm10_chromHMM_subfamily[which(mm10_chromHMM_subfamily$Mouse_sample_chromHMM == x[3]),],by="subfamily"))
hg19_mm10_chromHMM_subfamily = ldply(hg19_mm10_chromHMM_subfamily)
hg19_mm10_chromHMM_subfamily$Human_state_chromHMM = factor(hg19_mm10_chromHMM_subfamily$Human_state_chromHMM,chromHMM_states)
hg19_mm10_chromHMM_subfamily$Mouse_state_chromHMM = factor(hg19_mm10_chromHMM_subfamily$Mouse_state_chromHMM,mouse_chromHMM_states)
### Add sample categories
hg19_mm10_chromHMM_subfamily = merge(hg19_mm10_chromHMM_subfamily,human_mouse_samples[,c("Human_sample","Mouse_sample_chromHMM",sample_categories)],by=c("Human_sample","Mouse_sample_chromHMM"))

rm(list=c("hg19_chromHMM_subfamily","mm10_chromHMM_subfamily"))


# Proportion of each shared subfamily in an active regulatory chromHMM state in each species in each human/mouse sample pair
## Proportion of each hg19 subfamily in an active regulatory chromHMM state per sample, excluding those with 30 or fewer members
hg19_chromHMM_subfamily_active = read.table("Mouse/chromHMM/hg19_chromHMM_subfamily_active.txt",sep='\t',col.names=c("subfamily","Human_sample","Members"))
### Add missing subfamily x sample combinations
hg19_chromHMM_subfamily_active = merge(hg19_chromHMM_subfamily_active,expand.grid(subfamily=as.vector(rmsk_TE_subfamily$subfamily),
                                                                                  Human_sample=as.vector(human_mouse_samples$Human_sample)),
                                       by=c("subfamily","Human_sample"),all.y=TRUE)
hg19_chromHMM_subfamily_active[is.na(hg19_chromHMM_subfamily_active)] = 0
### Add number of subfamily members
hg19_chromHMM_subfamily_active = merge(hg19_chromHMM_subfamily_active,rmsk_TE_subfamily[,c("subfamily","Count")],by="subfamily")
### Calculate proprotion of subfamily members in state per sample
hg19_chromHMM_subfamily_active$Percent = hg19_chromHMM_subfamily_active$Members/hg19_chromHMM_subfamily_active$Count
### Restrict to subfamilies with > 30 members 
hg19_chromHMM_subfamily_active = hg19_chromHMM_subfamily_active[which(hg19_chromHMM_subfamily_active$Count > THRESHOLD_IK_MEMBER),]

## Proportion of each mm10 subfamily in an active regulatory chromHMM state per sample, excluding those with 30 or fewer members
mm10_chromHMM_subfamily_active = read.table("Mouse/chromHMM/mm10_chromHMM_subfamily_active.txt",sep='\t',col.names=c("subfamily","Mouse_sample_chromHMM","Members"))
### Add missing subfamily x sample combinations
mm10_chromHMM_subfamily_active = merge(mm10_chromHMM_subfamily_active,expand.grid(subfamily=as.vector(unique(mm10_rmsk_TE$subfamily)),
                                                                                  Mouse_sample_chromHMM=as.vector(human_mouse_samples$Mouse_sample_chromHMM)),
                                by=c("subfamily","Mouse_sample_chromHMM"),all.y=TRUE)
mm10_chromHMM_subfamily_active[is.na(mm10_chromHMM_subfamily_active)] = 0
### Add number of subfamily members
mm10_chromHMM_subfamily_active = merge(mm10_chromHMM_subfamily_active,mm10_subfamily_count[,c("subfamily","Count")],by="subfamily")
### Calculate proprotion of subfamily members in state per sample
mm10_chromHMM_subfamily_active$Percent = mm10_chromHMM_subfamily_active$Members/mm10_chromHMM_subfamily_active$Count
### Restrict to subfamilies with > 30 members 
mm10_chromHMM_subfamily_active = mm10_chromHMM_subfamily_active[which(mm10_chromHMM_subfamily_active$Count > THRESHOLD_IK_MEMBER),]

## Combine, restricting to human/mouse sample pairs, shared subfamilies only
hg19_mm10_chromHMM_subfamily_active = apply(human_mouse_samples,1,function(x) merge(hg19_chromHMM_subfamily_active[which(hg19_chromHMM_subfamily_active$Human_sample == x[1]),],
                                                                                                                  mm10_chromHMM_subfamily_active[which(mm10_chromHMM_subfamily_active$Mouse_sample_chromHMM == x[3]),],by="subfamily"))
hg19_mm10_chromHMM_subfamily_active = ldply(hg19_mm10_chromHMM_subfamily_active)
### Add sample categories
hg19_mm10_chromHMM_subfamily_active = merge(hg19_mm10_chromHMM_subfamily_active,human_mouse_samples[,c("Human_sample","Mouse_sample_chromHMM",sample_categories)],by=c("Human_sample","Mouse_sample_chromHMM"))
## Calculate the ratio of members in an active regulatory state in human versus mouse for each subfamily
hg19_mm10_chromHMM_subfamily_active$Ratio = hg19_mm10_chromHMM_subfamily_active$Percent.x/hg19_mm10_chromHMM_subfamily_active$Percent.y

rm(list=c("hg19_chromHMM_subfamily_active","mm10_chromHMM_subfamily_active","mm10_subfamily_count"))