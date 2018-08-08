# Subfamily analysis
#source("R_scripts/WGBS_subfamily_enrichment_TE.R")

# WGBS
## Proportion of mm10 subfamily hypomethylated
mm10_rmsk_TE_WGBS_subfamily_hypo = aggregate(data=mm10_rmsk_TE_WGBS[,c("subfamily",as.vector(na.omit(human_mouse_samples)$Mouse_sample_WGBS))],.~subfamily,function(x) sum(na.omit(x) < 0.3)/length(x),na.action=na.pass)

## Proportion of hg19 subfamily hypomethylated
hg19_rmsk_TE_WGBS_subfamily_hypo = TE_meth_subfamily[["Hypomethylated"]][,c("subfamily","family","class_update",as.vector(na.omit(human_mouse_samples)$Human_sample))]
hg19_rmsk_TE_WGBS_subfamily_hypo[,4:10] = t(apply(hg19_rmsk_TE_WGBS_subfamily_hypo,1,function(x) as.numeric(x[4:10])/rmsk_TE_subfamily[match(x[1],rmsk_TE_subfamily$subfamily),]$Count_CpG))

## Proportion of mm10 and hg19 subfamily hypomethylated
hg19_mm10_TE_WGBS_subfamily_hypo = merge(hg19_rmsk_TE_WGBS_subfamily_hypo,mm10_rmsk_TE_WGBS_subfamily_hypo,by="subfamily")
rm(list=c("hg19_rmsk_TE_WGBS_subfamily_hypo","mm10_rmsk_TE_WGBS_subfamily_hypo"))

# Filter out small subfamilies
mm10_subfamily_count = ddply(mm10_rmsk_TE,.(subfamily),summarise,TEs=length(subfamily))
hg19_mm10_TE_WGBS_subfamily_hypo = hg19_mm10_TE_WGBS_subfamily_hypo[which(!(hg19_mm10_TE_WGBS_subfamily_hypo$subfamily %in% 
                                                                                            union(as.vector(rmsk_TE_subfamily[which(rmsk_TE_subfamily$Count <= THRESHOLD_IK_MEMBER),]$subfamily),
                                                                                                  as.vector(mm10_subfamily_count[which(mm10_subfamily_count$TEs <= THRESHOLD_IK_MEMBER),]$subfamily)))),]
## Paired samples only
hg19_mm10_TE_WGBS_subfamily_hypo_paired = melt(hg19_mm10_TE_WGBS_subfamily_hypo,id.vars=c("subfamily","family","class_update",as.vector(na.omit(human_mouse_samples)$Human_sample)))
colnames(hg19_mm10_TE_WGBS_subfamily_hypo_paired)[11:12] = c("Mouse_sample_WGBS","Mouse_hypo")
hg19_mm10_TE_WGBS_subfamily_hypo_paired = melt(hg19_mm10_TE_WGBS_subfamily_hypo_paired,id.vars=c("subfamily","family","class_update","Mouse_sample_WGBS","Mouse_hypo"))
colnames(hg19_mm10_TE_WGBS_subfamily_hypo_paired)[6:7]= c("Human_sample","Human_hypo")
hg19_mm10_TE_WGBS_subfamily_hypo_paired = merge(hg19_mm10_TE_WGBS_subfamily_hypo_paired,human_mouse_samples,by=c("Mouse_sample_WGBS","Human_sample"))

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

## mm10 subfamily proportion in state
mm10_chromHMM_subfamily = read.table("Mouse/chromHMM/mm10_chromHMM_subfamily.txt",sep = '\t')
colnames(mm10_chromHMM_subfamily) = c("subfamily","Mouse_state_chromHMM","Mouse_sample_chromHMM","Members")
mm10_chromHMM_subfamily = merge(mm10_chromHMM_subfamily,
                                expand.grid(subfamily=as.vector(unique(mm10_rmsk_TE$subfamily)),Mouse_state_chromHMM=mouse_chromHMM_states,Mouse_sample_chromHMM=as.vector(human_mouse_samples$Mouse_sample_chromHMM)),
                                by=c("subfamily","Mouse_state_chromHMM","Mouse_sample_chromHMM"),all.y=TRUE)
mm10_chromHMM_subfamily[is.na(mm10_chromHMM_subfamily)] = 0
mm10_chromHMM_subfamily = merge(mm10_chromHMM_subfamily,ddply(mm10_rmsk_TE,.(subfamily),summarise,Count=length(subfamily)),by="subfamily",all.y=TRUE)
mm10_chromHMM_subfamily$Percent = mm10_chromHMM_subfamily$Members/mm10_chromHMM_subfamily$Count
mm10_chromHMM_subfamily = mm10_chromHMM_subfamily[which(mm10_chromHMM_subfamily$Count > THRESHOLD_IK_MEMBER),]

# Spearman correlation between all TE subfamilies for each sample pair, chromHMM
## Promoter
human_mouse_subfamily_TssA = merge(dcast(hg19_chromHMM_subfamily[which(hg19_chromHMM_subfamily$State == "1_TssA"),],subfamily~Sample,value.var="Percent"),
             dcast(mm10_chromHMM_subfamily[which(mm10_chromHMM_subfamily$Mouse_state_chromHMM == "TssA"),],subfamily~Mouse_sample_chromHMM,value.var="Percent"),by="subfamily")
human_mouse_spearman_subfamily_TssA = apply(human_mouse_subfamily_TssA[,c(as.vector(human_mouse_samples$Mouse_sample_chromHMM),as.vector(human_mouse_samples$Human_sample))],2,function(x) 
  apply(human_mouse_subfamily_TssA[,c(as.vector(human_mouse_samples$Mouse_sample_chromHMM),as.vector(human_mouse_samples$Human_sample))],2,function(y) unlist(cor.test(x,y,method="spearman"))["estimate.rho"]))
human_mouse_spearman_subfamily_TssA = melt(as.matrix(human_mouse_spearman_subfamily_TssA))
colnames(human_mouse_spearman_subfamily_TssA) = c("Sample1","Sample2","Correlation")
human_mouse_spearman_subfamily_TssA = merge(human_mouse_spearman_subfamily_TssA,human_mouse_samples_expand)
human_mouse_spearman_subfamily_TssA$Correlation = as.numeric(as.character(human_mouse_spearman_subfamily_TssA$Correlation))

### Paired samples only
human_mouse_subfamily_TssA_paired = melt(human_mouse_subfamily_TssA,id.vars=c("subfamily",as.vector(human_mouse_samples$Human_sample)))
colnames(human_mouse_subfamily_TssA_paired)[14:15] = c("Mouse_sample_chromHMM","Mouse_proportion")
human_mouse_subfamily_TssA_paired = melt(human_mouse_subfamily_TssA_paired,id.vars=c("subfamily","Mouse_sample_chromHMM","Mouse_proportion"))
colnames(human_mouse_subfamily_TssA_paired)[4:5]= c("Human_sample","Human_proportion")
human_mouse_subfamily_TssA_paired = merge(human_mouse_subfamily_TssA_paired,human_mouse_samples,by=c("Mouse_sample_chromHMM","Human_sample"))
human_mouse_subfamily_TssA_paired$State = rep("1_TssA",dim(human_mouse_subfamily_TssA_paired)[1])

## Enhancer
human_mouse_subfamily_Enh = merge(dcast(hg19_chromHMM_subfamily[which(hg19_chromHMM_subfamily$State == "7_Enh"),],subfamily~Sample,value.var="Percent"),
             dcast(mm10_chromHMM_subfamily[which(mm10_chromHMM_subfamily$Mouse_state_chromHMM == "Enh"),],subfamily~Mouse_sample_chromHMM,value.var="Percent"),by="subfamily")
human_mouse_spearman_subfamily_Enh = apply(human_mouse_subfamily_Enh[,c(as.vector(human_mouse_samples$Mouse_sample_chromHMM),as.vector(human_mouse_samples$Human_sample))],2,function(x) 
  apply(human_mouse_subfamily_Enh[,c(as.vector(human_mouse_samples$Mouse_sample_chromHMM),as.vector(human_mouse_samples$Human_sample))],2,function(y) unlist(cor.test(x,y,method="spearman"))["estimate.rho"]))
human_mouse_spearman_subfamily_Enh = melt(as.matrix(human_mouse_spearman_subfamily_Enh))
colnames(human_mouse_spearman_subfamily_Enh) = c("Sample1","Sample2","Correlation")
human_mouse_spearman_subfamily_Enh = merge(human_mouse_spearman_subfamily_Enh,human_mouse_samples_expand)
human_mouse_spearman_subfamily_Enh$Correlation = as.numeric(as.character(human_mouse_spearman_subfamily_Enh$Correlation))

### Paired samples only
human_mouse_subfamily_Enh_paired = melt(human_mouse_subfamily_Enh,id.vars=c("subfamily",as.vector(human_mouse_samples$Human_sample)))
colnames(human_mouse_subfamily_Enh_paired)[14:15] = c("Mouse_sample_chromHMM","Mouse_proportion")
human_mouse_subfamily_Enh_paired = melt(human_mouse_subfamily_Enh_paired,id.vars=c("subfamily","Mouse_sample_chromHMM","Mouse_proportion"))
colnames(human_mouse_subfamily_Enh_paired)[4:5]= c("Human_sample","Human_proportion")
human_mouse_subfamily_Enh_paired = merge(human_mouse_subfamily_Enh_paired,human_mouse_samples,by=c("Mouse_sample_chromHMM","Human_sample"))
human_mouse_subfamily_Enh_paired$State = rep("7_Enh",dim(human_mouse_subfamily_Enh_paired)[1])