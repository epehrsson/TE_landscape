# Comparison of human-mouse TE ortholog methylation state
# See 5/11/2017, 5/12/2017, 5/22/2017, 5/29/2017, 5/30/2017, 7/21/2017, 7/26/2017, 7/31/2017, 8/1/2017

# WGBS
## hg19 TEs with WGBS average methylation

## mm10 TEs with WGBS average methylation
mm10_rmsk_TE_WGBS = read.table("Mouse/WGBS/mm10_rmsk_TE_WGBS_avg.txt",sep='\t',header=TRUE)
mm10_rmsk_TE_WGBS[,8:16] = mm10_rmsk_TE_WGBS[,8:16]/100
mm10_rmsk_TE_CpG_count = read.table("Mouse/WGBS/mm10_rmsk_TE_CpG_count.txt",sep='\t')
colnames(mm10_rmsk_TE_CpG_count) = c(TE_coordinates[c(1:4,6,5,7)],"CpGs")
mm10_rmsk_TE_CpG_count$CpGs = mm10_rmsk_TE_CpG_count$CpGs/2
mm10_rmsk_TE_WGBS = merge(mm10_rmsk_TE_WGBS,mm10_rmsk_TE_CpG_count,by=TE_coordinates)
rm(mm10_rmsk_TE_CpG_count)

## Orthologous TEs with methylation level (only those with CpGs in both)
human_mouse_orthologs_WGBS = merge(human_mouse_orthologs_mm10[,c(hg19_coordinates,mm10_coordinates,"class_update","JC_distance")],
                                   TE_meth_average[,c(TE_coordinates,as.vector(na.omit(human_mouse_samples)$Human_sample))],by.y=TE_coordinates[c(1:4,6,5,7)],by.x=hg19_coordinates)
human_mouse_orthologs_WGBS = merge(human_mouse_orthologs_WGBS,mm10_rmsk_TE_WGBS[,c(TE_coordinates,as.vector(na.omit(human_mouse_samples)$Mouse_sample_WGBS))],by.y=TE_coordinates[c(1:4,6,5,7)],by.x=mm10_coordinates)

## Paired samples only
human_mouse_orthologs_WGBS_paired = melt(human_mouse_orthologs_WGBS,id.vars=c(hg19_coordinates,mm10_coordinates,"class_update","JC_distance",as.vector(na.omit(human_mouse_samples)$Human)))
colnames(human_mouse_orthologs_WGBS_paired)[24:25] = c("Mouse_sample_WGBS","Mouse_methylation")
human_mouse_orthologs_WGBS_paired = melt(human_mouse_orthologs_WGBS_paired,id.vars=c(hg19_coordinates,mm10_coordinates,"class_update","JC_distance","Mouse_sample_WGBS","Mouse_methylation"))
colnames(human_mouse_orthologs_WGBS_paired)[19:20] = c("Human_sample","Human_methylation")
human_mouse_orthologs_WGBS_paired = merge(human_mouse_orthologs_WGBS_paired,na.omit(human_mouse_samples)[,c("Mouse_sample_WGBS","Human_sample")],by=c("Mouse_sample_WGBS","Human_sample"))

## Adding methylation state
human_mouse_orthologs_WGBS_paired$Human_state_WGBS = factor(ifelse(is.na(human_mouse_orthologs_WGBS_paired$Human_methylation),"Missing",
                                                       ifelse(human_mouse_orthologs_WGBS_paired$Human_methylation > 0.7,"Hypermethylated",
                                                       ifelse(human_mouse_orthologs_WGBS_paired$Human_methylation < 0.3,"Hypomethylated",
                                                              ifelse(human_mouse_orthologs_WGBS_paired$Human_methylation <= 0.7 & human_mouse_orthologs_WGBS_paired$Human_methylation >= 0.3,"Intermediate",NA)))),
                                                       levels=meth_states)
human_mouse_orthologs_WGBS_paired$Mouse_state_WGBS = factor(ifelse(is.na(human_mouse_orthologs_WGBS_paired$Mouse_methylation),"Missing",
                                                       ifelse(human_mouse_orthologs_WGBS_paired$Mouse_methylation > 0.7,"Hypermethylated",
                                                       ifelse(human_mouse_orthologs_WGBS_paired$Mouse_methylation < 0.3,"Hypomethylated",
                                                              ifelse(human_mouse_orthologs_WGBS_paired$Mouse_methylation <= 0.7 & human_mouse_orthologs_WGBS_paired$Mouse_methylation >= 0.3,"Intermediate",NA)))),
                                                       levels=meth_states)

# Spearman correlation between all TEs for each sample pair
human_mouse_spearman = apply(human_mouse_orthologs_WGBS[,c(as.vector(na.omit(human_mouse_samples)$Mouse_sample_WGBS),as.vector(na.omit(human_mouse_samples)$Human_sample))],2,function(x)
  apply(human_mouse_orthologs_WGBS[,c(as.vector(na.omit(human_mouse_samples)$Mouse_sample_WGBS),as.vector(na.omit(human_mouse_samples)$Human_sample))],2,function(y) unlist(cor.test(x,y,method="spearman"))["estimate.rho"]))
human_mouse_spearman = melt(as.matrix(human_mouse_spearman))
colnames(human_mouse_spearman) = c("Sample1","Sample2","Correlation")
human_mouse_spearman = merge(human_mouse_spearman,human_mouse_samples_WGBS)
human_mouse_spearman$Correlation = as.numeric(as.character(human_mouse_spearman$Correlation))