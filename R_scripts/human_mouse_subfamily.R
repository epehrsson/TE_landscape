# Subfamily analysis
#source("R_scripts/WGBS_subfamily_enrichment_TE.R")

# WGBS
## Proportion of mm10 subfamily hypomethylated
mm10_rmsk_TE_WGBS_subfamily_hypo = aggregate(data=mm10_rmsk_TE_WGBS[,c("subfamily",as.vector(na.omit(human_mouse_samples)$Mouse_sample_WGBS))],.~subfamily,function(x) sum(na.omit(x) < 0.3)/length(x),na.action=na.pass)

## Proportion of hg19 subfamily hypomethylated
hg19_rmsk_TE_WGBS_subfamily_hypo = merge(TE_meth_subfamily[which(TE_meth_subfamily$State == "Hypomethylated" & TE_meth_subfamily$Sample %in% as.vector(na.omit(human_mouse_samples)$Human_sample)),
                                                     c("subfamily","family","class_update","Sample","Members")],rmsk_TE_subfamily[,c("subfamily","Count_CpGs")],by="subfamily")
hg19_rmsk_TE_WGBS_subfamily_hypo$Percent = hg19_rmsk_TE_WGBS_subfamily_hypo$Members/hg19_rmsk_TE_WGBS_subfamily_hypo$Count_CpGs
hg19_rmsk_TE_WGBS_subfamily_hypo = dcast(hg19_rmsk_TE_WGBS_subfamily_hypo,subfamily+family+class_update~Sample,value.var="Percent")

## Proportion of mm10 and hg19 subfamily hypomethylated
hg19_mm10_TE_WGBS_subfamily_hypo = merge(hg19_rmsk_TE_WGBS_subfamily_hypo,mm10_rmsk_TE_WGBS_subfamily_hypo,by="subfamily")
rm(list=c("hg19_rmsk_TE_WGBS_subfamily_hypo","mm10_rmsk_TE_WGBS_subfamily_hypo"))

# Filter out small subfamilies (TEs with CpGs)
mm10_subfamily_count = merge(ddply(mm10_rmsk_TE,.(subfamily),summarise,Count=length(subfamily)),ddply(mm10_rmsk_TE_WGBS,.(subfamily),summarise,Count_CpGs=length(subfamily)),by="subfamily",all.x=TRUE)
mm10_subfamily_count[is.na(mm10_subfamily_count)] = 0
hg19_mm10_subfamily_count = merge(rmsk_TE_subfamily[,c("subfamily","family","class_update","Count","Count_CpGs")],mm10_subfamily_count,by="subfamily")

hg19_mm10_TE_WGBS_subfamily_hypo = hg19_mm10_TE_WGBS_subfamily_hypo[which(hg19_mm10_TE_WGBS_subfamily_hypo$subfamily %in% as.vector(hg19_mm10_subfamily_count[which(hg19_mm10_subfamily_count$Count_CpGs.x > THRESHOLD_IK_MEMBER & 
                                                                                                                                                                    hg19_mm10_subfamily_count$Count_CpGs.y > THRESHOLD_IK_MEMBER),]$subfamily)),]

## Paired samples only
hg19_mm10_TE_WGBS_subfamily_hypo_paired = melt(hg19_mm10_TE_WGBS_subfamily_hypo,id.vars=c("subfamily","family","class_update",as.vector(na.omit(human_mouse_samples)$Human_sample)))
colnames(hg19_mm10_TE_WGBS_subfamily_hypo_paired)[11:12] = c("Mouse_sample_WGBS","Mouse_hypo")
hg19_mm10_TE_WGBS_subfamily_hypo_paired = melt(hg19_mm10_TE_WGBS_subfamily_hypo_paired,id.vars=c("subfamily","family","class_update","Mouse_sample_WGBS","Mouse_hypo"))
colnames(hg19_mm10_TE_WGBS_subfamily_hypo_paired)[6:7]= c("Human_sample","Human_hypo")
hg19_mm10_TE_WGBS_subfamily_hypo_paired = merge(hg19_mm10_TE_WGBS_subfamily_hypo_paired,human_mouse_samples,by=c("Mouse_sample_WGBS","Human_sample"))
hg19_mm10_TE_WGBS_subfamily_hypo_paired$Ratio = hg19_mm10_TE_WGBS_subfamily_hypo_paired$Human_hypo/hg19_mm10_TE_WGBS_subfamily_hypo_paired$Mouse_hypo

# Spearman correlation between all TE subfamilies for each sample pair, WGBS
human_mouse_spearman_subfamily = apply(hg19_mm10_TE_WGBS_subfamily_hypo[,c(as.vector(na.omit(human_mouse_samples)$Mouse_sample_WGBS),as.vector(na.omit(human_mouse_samples)$Human_sample))],2,function(x) 
  apply(hg19_mm10_TE_WGBS_subfamily_hypo[,c(as.vector(na.omit(human_mouse_samples)$Mouse_sample_WGBS),as.vector(na.omit(human_mouse_samples)$Human_sample))],2,function(y) unlist(cor.test(x,y,method="spearman"))["estimate.rho"]))
human_mouse_spearman_subfamily = melt(as.matrix(human_mouse_spearman_subfamily))
colnames(human_mouse_spearman_subfamily) = c("Sample1","Sample2","Correlation")
human_mouse_spearman_subfamily = merge(human_mouse_spearman_subfamily,na.omit(human_mouse_samples_WGBS))
human_mouse_spearman_subfamily$Correlation = as.numeric(as.character(human_mouse_spearman_subfamily$Correlation))

# chromHMM
## hg19 subfamily proportion in state (removes small subfamilies)
hg19_chromHMM_subfamily = subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% chromHMM_states & subfamily_state_sample_combined$Sample %in% as.vector(human_mouse_samples$Human_sample) & 
                                                                  subfamily_state_sample_combined$Count > THRESHOLD_IK_MEMBER),
                                                          c("subfamily","family","class_update","Sample","State","Members","Count","Percent")]
colnames(hg19_chromHMM_subfamily)[4:5] = c("Human_sample","Human_state_chromHMM")

## mm10 subfamily proportion in state (removes small subfamilies)
mm10_chromHMM_subfamily = read.table("Mouse/chromHMM/mm10_chromHMM_subfamily.txt",sep = '\t')
colnames(mm10_chromHMM_subfamily) = c("subfamily","Mouse_state_chromHMM","Mouse_sample_chromHMM","Members")
mm10_chromHMM_subfamily = merge(mm10_chromHMM_subfamily,
                                expand.grid(subfamily=as.vector(unique(mm10_rmsk_TE$subfamily)),Mouse_state_chromHMM=mouse_chromHMM_states,Mouse_sample_chromHMM=as.vector(human_mouse_samples$Mouse_sample_chromHMM)),
                                by=c("subfamily","Mouse_state_chromHMM","Mouse_sample_chromHMM"),all.y=TRUE)
mm10_chromHMM_subfamily[is.na(mm10_chromHMM_subfamily)] = 0
mm10_chromHMM_subfamily = merge(mm10_chromHMM_subfamily,mm10_subfamily_count[,c("subfamily","Count")],by="subfamily")
mm10_chromHMM_subfamily$Percent = mm10_chromHMM_subfamily$Members/mm10_chromHMM_subfamily$Count
mm10_chromHMM_subfamily = mm10_chromHMM_subfamily[which(mm10_chromHMM_subfamily$Count > THRESHOLD_IK_MEMBER),]

## Proportion of mm10 and hg19 subfamily in each chromHMM state, paired samples only
hg19_mm10_chromHMM_subfamily = apply(human_mouse_samples,1,function(x) merge(hg19_chromHMM_subfamily[which(hg19_chromHMM_subfamily$Human_sample == x[1]),],
                                                                             mm10_chromHMM_subfamily[which(mm10_chromHMM_subfamily$Mouse_sample_chromHMM == x[3]),],by="subfamily"))
hg19_mm10_chromHMM_subfamily = ldply(hg19_mm10_chromHMM_subfamily)
hg19_mm10_chromHMM_subfamily$Human_state_chromHMM = factor(hg19_mm10_chromHMM_subfamily$Human_state_chromHMM,chromHMM_states)
hg19_mm10_chromHMM_subfamily$Mouse_state_chromHMM = factor(hg19_mm10_chromHMM_subfamily$Mouse_state_chromHMM,mouse_chromHMM_states)
hg19_mm10_chromHMM_subfamily = merge(hg19_mm10_chromHMM_subfamily,human_mouse_samples[,c("Human_sample","Mouse_sample_chromHMM",sample_categories)],by=c("Human_sample","Mouse_sample_chromHMM"))

## Spearman correlation between all TE subfamilies for each sample pair, promoter
human_mouse_subfamily_TssA = merge(dcast(hg19_chromHMM_subfamily[which(hg19_chromHMM_subfamily$Human_state_chromHMM == "1_TssA"),],subfamily~Human_sample,value.var="Percent"),
             dcast(mm10_chromHMM_subfamily[which(mm10_chromHMM_subfamily$Mouse_state_chromHMM == "TssA"),],subfamily~Mouse_sample_chromHMM,value.var="Percent"),by="subfamily")
human_mouse_spearman_subfamily_TssA = apply(human_mouse_subfamily_TssA[,c(as.vector(human_mouse_samples$Mouse_sample_chromHMM),as.vector(human_mouse_samples$Human_sample))],2,function(x) 
  apply(human_mouse_subfamily_TssA[,c(as.vector(human_mouse_samples$Mouse_sample_chromHMM),as.vector(human_mouse_samples$Human_sample))],2,function(y) unlist(cor.test(x,y,method="spearman"))["estimate.rho"]))
human_mouse_spearman_subfamily_TssA = melt(as.matrix(human_mouse_spearman_subfamily_TssA))
colnames(human_mouse_spearman_subfamily_TssA) = c("Sample1","Sample2","Correlation")
human_mouse_spearman_subfamily_TssA = merge(human_mouse_spearman_subfamily_TssA,human_mouse_samples_expand)
human_mouse_spearman_subfamily_TssA$Correlation = as.numeric(as.character(human_mouse_spearman_subfamily_TssA$Correlation))

# Proportion of each subfamily in active regulatory states, paired samples only
## hg19 subfamily proportion
hg19_chromHMM_subfamily_active = read.table("Mouse/chromHMM/hg19_chromHMM_subfamily_active.txt",sep='\t')
colnames(hg19_chromHMM_subfamily_active) = c("subfamily","Human_sample","Members")
hg19_chromHMM_subfamily_active = merge(hg19_chromHMM_subfamily_active,expand.grid(subfamily=as.vector(rmsk_TE_subfamily$subfamily),Human_sample=as.vector(human_mouse_samples$Human_sample)),
                                       by=c("subfamily","Human_sample"),all.y=TRUE)
hg19_chromHMM_subfamily_active[is.na(hg19_chromHMM_subfamily_active)] = 0
hg19_chromHMM_subfamily_active = merge(hg19_chromHMM_subfamily_active,rmsk_TE_subfamily[,c("subfamily","Count")],by="subfamily")
hg19_chromHMM_subfamily_active$Percent = hg19_chromHMM_subfamily_active$Members/hg19_chromHMM_subfamily_active$Count
hg19_chromHMM_subfamily_active = hg19_chromHMM_subfamily_active[which(hg19_chromHMM_subfamily_active$Count > THRESHOLD_IK_MEMBER),]

## mm10 subfamily proportion
mm10_chromHMM_subfamily_active = read.table("Mouse/chromHMM/mm10_chromHMM_subfamily_active.txt",sep='\t')
colnames(mm10_chromHMM_subfamily_active) = c("subfamily","Mouse_sample_chromHMM","Members")
mm10_chromHMM_subfamily_active = merge(mm10_chromHMM_subfamily_active,expand.grid(subfamily=as.vector(unique(mm10_rmsk_TE$subfamily)),Mouse_sample_chromHMM=as.vector(human_mouse_samples$Mouse_sample_chromHMM)),
                                by=c("subfamily","Mouse_sample_chromHMM"),all.y=TRUE)
mm10_chromHMM_subfamily_active[is.na(mm10_chromHMM_subfamily_active)] = 0
mm10_chromHMM_subfamily_active = merge(mm10_chromHMM_subfamily_active,mm10_subfamily_count[,c("subfamily","Count")],by="subfamily")
mm10_chromHMM_subfamily_active$Percent = mm10_chromHMM_subfamily_active$Members/mm10_chromHMM_subfamily_active$Count
mm10_chromHMM_subfamily_active = mm10_chromHMM_subfamily_active[which(mm10_chromHMM_subfamily_active$Count > THRESHOLD_IK_MEMBER),]

## Paired samples
hg19_mm10_chromHMM_subfamily_active = apply(human_mouse_samples,1,function(x) merge(hg19_chromHMM_subfamily_active[which(hg19_chromHMM_subfamily_active$Human_sample == x[1]),],
                                                                                                                  mm10_chromHMM_subfamily_active[which(mm10_chromHMM_subfamily_active$Mouse_sample_chromHMM == x[3]),],by="subfamily"))
hg19_mm10_chromHMM_subfamily_active = ldply(hg19_mm10_chromHMM_subfamily_active)
hg19_mm10_chromHMM_subfamily_active = merge(hg19_mm10_chromHMM_subfamily_active,human_mouse_samples[,c("Human_sample","Mouse_sample_chromHMM",sample_categories)],by=c("Human_sample","Mouse_sample_chromHMM"))
hg19_mm10_chromHMM_subfamily_active$Ratio = hg19_mm10_chromHMM_subfamily_active$Percent.x/hg19_mm10_chromHMM_subfamily_active$Percent.y

## Spearman correlation between all TE subfamilies for each sample pair, active regulatory
human_mouse_subfamily_active = merge(dcast(hg19_chromHMM_subfamily_active,subfamily~Human_sample,value.var="Percent"),
             dcast(mm10_chromHMM_subfamily_active,subfamily~Mouse_sample_chromHMM,value.var="Percent"),by="subfamily")
human_mouse_spearman_subfamily_active = apply(human_mouse_subfamily_active[,c(as.vector(human_mouse_samples$Mouse_sample_chromHMM),as.vector(human_mouse_samples$Human_sample))],2,function(x) 
  apply(human_mouse_subfamily_active[,c(as.vector(human_mouse_samples$Mouse_sample_chromHMM),as.vector(human_mouse_samples$Human_sample))],2,function(y) unlist(cor.test(x,y,method="spearman"))["estimate.rho"]))
human_mouse_spearman_subfamily_active = melt(as.matrix(human_mouse_spearman_subfamily_active))
colnames(human_mouse_spearman_subfamily_active) = c("Sample1","Sample2","Correlation")
human_mouse_spearman_subfamily_active = merge(human_mouse_spearman_subfamily_active,human_mouse_samples_expand)
human_mouse_spearman_subfamily_active$Correlation = as.numeric(as.character(human_mouse_spearman_subfamily_active$Correlation))