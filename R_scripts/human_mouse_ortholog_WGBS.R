# Generates dataframes of the average methylation level and methylation state for orthologous TEs,
# including mm10 TEs and hg19/mm10 ortholog pairs

## mm10_rmsk_TE_WGBS - Dataframe of average methylation for each mm10 TE that overlaps a CpG, in 9 mouse samples
## human_mouse_orthologs_WGBS - Dataframe of orthologous TE (hg19/mm10) average methylation in human and mouse samples, 
## including only TEs that overlap a CpG in each species
## human_mouse_orthologs_WGBS_paired - Dataframe of orthologous TE average methylation in human and mouse samples, restricted to paired samples


# Dataframe of average methylation for each mm10 TE in 9 mouse samples
mm10_rmsk_TE_WGBS = read.table("Mouse/WGBS/mm10_rmsk_TE_WGBS_avg.txt",sep='\t',header=TRUE)
## Convert methylation from percent to proportion
mm10_rmsk_TE_WGBS[,8:16] = mm10_rmsk_TE_WGBS[,8:16]/100

## Load number of CpGs per mm10 TE
mm10_rmsk_TE_CpG_count = read.table("Mouse/WGBS/mm10_rmsk_TE_CpG_count.txt",sep='\t',col.names=c(TE_coordinates[c(1:4,6,5,7)],"CpGs"))
mm10_rmsk_TE_CpG_count$CpGs = mm10_rmsk_TE_CpG_count$CpGs/2

## Restrict TEs to only those overlapping any CpG
mm10_rmsk_TE_WGBS = merge(mm10_rmsk_TE_WGBS,mm10_rmsk_TE_CpG_count,by=TE_coordinates)
rm(mm10_rmsk_TE_CpG_count)


# Dataframe of orthologous TEs (hg19/mm10) with average methylation in human and mouse samples, including only TEs that overlap a CpG in each species
## Add average methylation for each hg19 TEs with an mm10 ortholog in human samples with a corresponding mouse ENCODE WGBS sample
human_mouse_orthologs_WGBS = merge(human_mouse_orthologs_mm10[,c(hg19_coordinates,mm10_coordinates,"class_update","JC_distance")],
                                   TE_meth_average[,c(TE_coordinates,as.vector(na.omit(human_mouse_samples)$Human_sample))],
                                   by.y=TE_coordinates[c(1:4,6,5,7)],by.x=hg19_coordinates)
## Add average methylation for each mm10 orthologous TE in mouse ENCODE WGBS samples
human_mouse_orthologs_WGBS = merge(human_mouse_orthologs_WGBS,mm10_rmsk_TE_WGBS[,c(TE_coordinates,as.vector(na.omit(human_mouse_samples)$Mouse_sample_WGBS))],
                                   by.y=TE_coordinates[c(1:4,6,5,7)],by.x=mm10_coordinates)


# Restrict dataframe of orthologous TE average methylation in human and mouse samples to paired samples (same tissue) and add state
## Reformat dataframe so that human/mouse samples and methylation are in long format
human_mouse_orthologs_WGBS_paired = melt(human_mouse_orthologs_WGBS,
                                         id.vars=c(hg19_coordinates,mm10_coordinates,"class_update","JC_distance",as.vector(na.omit(human_mouse_samples)$Human)),
                                         variable.name="Mouse_sample_WGBS",value.name="Mouse_methylation")
human_mouse_orthologs_WGBS_paired = melt(human_mouse_orthologs_WGBS_paired,
                                         id.vars=c(hg19_coordinates,mm10_coordinates,"class_update","JC_distance","Mouse_sample_WGBS","Mouse_methylation"),
                                         variable.name="Human_sample",value.name="Human_methylation")
## Restrict dataframe to only corresponding human/mouse samples
human_mouse_orthologs_WGBS_paired = merge(human_mouse_orthologs_WGBS_paired,na.omit(human_mouse_samples)[,c("Mouse_sample_WGBS","Human_sample")],
                                          by=c("Mouse_sample_WGBS","Human_sample"))

## Convert average methylation to methylation state
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