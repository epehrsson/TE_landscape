# Subfamily analysis
# Proportion of mm10 subfamily hypomethylated
mm10_rmsk_TE_WGBS_subfamily_hypo = aggregate(data=mm10_rmsk_TE_WGBS[,c("subfamily",as.vector(unique(human_mouse_samples$Mouse)))],.~subfamily,function(x) sum(na.omit(x) < 0.3)/length(x),na.action=na.pass)

# Proportion of hg19 subfamily hypomethylated
hg19_rmsk_TE_WGBS_subfamily_hypo = TE_meth_subfamily[["Hypomethylated"]][,c("subfamily","family","class_update",as.vector(human_mouse_samples$Human))]
hg19_rmsk_TE_WGBS_subfamily_hypo[,4:16] = t(apply(hg19_rmsk_TE_WGBS_subfamily_hypo,1,function(x) as.numeric(x[4:16])/rmsk_TE_subfamily[match(x[1],rmsk_TE_subfamily$subfamily),]$Count_CpG))

# Proportion of mm10 and hg19 subfamily hypomethylated
hg19_mm10_TE_WGBS_subfamily_hypo = merge(hg19_rmsk_TE_WGBS_subfamily_hypo,mm10_rmsk_TE_WGBS_subfamily_hypo,by="subfamily")
rm(list=c("hg19_rmsk_TE_WGBS_subfamily_hypo","mm10_rmsk_TE_WGBS_subfamily_hypo"))

# Paired samples only
hg19_mm10_TE_WGBS_subfamily_hypo_paired = melt(hg19_mm10_TE_WGBS_subfamily_hypo,id.vars=c("subfamily","family","class_update",as.vector(human_mouse_samples$Human)))
colnames(hg19_mm10_TE_WGBS_subfamily_hypo_paired)[17:18] = c("Mouse","Mouse_proportion")
hg19_mm10_TE_WGBS_subfamily_hypo_paired = melt(hg19_mm10_TE_WGBS_subfamily_hypo_paired,id.vars=c("subfamily","family","class_update","Mouse","Mouse_proportion"))
colnames(hg19_mm10_TE_WGBS_subfamily_hypo_paired)[6:7]= c("Human","Human_proportion")
hg19_mm10_TE_WGBS_subfamily_hypo_paired = merge(hg19_mm10_TE_WGBS_subfamily_hypo_paired,human_mouse_samples,by=c("Mouse","Human"))

# Spearman correlation between all TE subfamilies for each sample pair
human_mouse_spearman_subfamily = apply(hg19_mm10_TE_WGBS_subfamily_hypo[,c(unique(as.vector(human_mouse_samples$Mouse)),as.vector(human_mouse_samples$Human))],2,function(x) 
  apply(hg19_mm10_TE_WGBS_subfamily_hypo[,c(unique(as.vector(human_mouse_samples$Mouse)),as.vector(human_mouse_samples$Human))],2,function(y) unlist(cor.test(x,y,method="spearman"))["estimate.rho"]))
human_mouse_spearman_subfamily = melt(as.matrix(human_mouse_spearman_subfamily))
colnames(human_mouse_spearman_subfamily) = c("Sample1","Sample2","Correlation")
human_mouse_spearman_subfamily = merge(human_mouse_spearman_subfamily,human_mouse_samples_expand)
human_mouse_spearman_subfamily$Correlation = as.numeric(as.character(human_mouse_spearman_subfamily$Correlation))