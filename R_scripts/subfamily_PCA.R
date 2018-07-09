# Perform PCA/tSNE at the subfamily level

set.seed(42)

# chromHMM 3D clustering - by sample
sample_chromHMM_3D = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% chromHMM_states),c("subfamily","State","Sample","Length_percent_jk")],
                                      subfamily+State~Sample,value.var = "Length_percent_jk")
rownames(sample_chromHMM_3D) = paste(sample_chromHMM_3D$subfamily,sample_chromHMM_3D$State,sep="_")
sample_chromHMM_3D = sample_chromHMM_3D[,3:129]
sample_chromHMM_3D = sample_chromHMM_3D[which(apply(sample_chromHMM_3D,1,var) != 0),]

# chromHMM 3D clustering - by subfamily
subfamily_chromHMM_3D = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% chromHMM_states),c("subfamily","State","Sample","Length_percent_jk")],
                           Sample+State~subfamily,value.var = "Length_percent_jk")
rownames(subfamily_chromHMM_3D) = paste(subfamily_chromHMM_3D$Sample,subfamily_chromHMM_3D$State,sep="_")
subfamily_chromHMM_3D = subfamily_chromHMM_3D[,3:970]
subfamily_chromHMM_3D = subfamily_chromHMM_3D[which(apply(subfamily_chromHMM_3D,1,var) != 0),]

# WGBS 3D clustering - by sample
sample_WGBS_3D = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% meth_states),c("subfamily","State","Sample","Length_percent_jk")],
                           subfamily+State~Sample,value.var = "Length_percent_jk")
rownames(sample_WGBS_3D) = paste(sample_WGBS_3D$subfamily,sample_WGBS_3D$State,sep="_")
sample_WGBS_3D = sample_WGBS_3D[,3:39]
sample_WGBS_3D = sample_WGBS_3D[which(apply(sample_WGBS_3D,1,var) != 0),]

# WGBS 3D clustering - by subfamily
subfamily_WGBS_3D = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% meth_states),c("subfamily","State","Sample","Length_percent_jk")],
                              Sample+State~subfamily,value.var = "Length_percent_jk")
rownames(subfamily_WGBS_3D) = paste(subfamily_WGBS_3D$Sample,subfamily_WGBS_3D$State,sep="_")
subfamily_WGBS_3D = subfamily_WGBS_3D[,3:967]
subfamily_WGBS_3D = subfamily_WGBS_3D[which(apply(subfamily_WGBS_3D,1,var) != 0),]

## DNase
subfamily_DNase_enrichment = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State == "DNase"),c("subfamily","Sample","Length_percent_jk")],
                                   subfamily~Sample,value.var = "Length_percent_jk")
rownames(subfamily_DNase_enrichment) = subfamily_DNase_enrichment$subfamily
subfamily_DNase_enrichment = subfamily_DNase_enrichment[,2:54]
subfamily_DNase_enrichment = subfamily_DNase_enrichment[which(apply(subfamily_DNase_enrichment,1,var) != 0),]

## H3K27ac
subfamily_H3K27ac_enrichment = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State == "H3K27ac"),c("subfamily","Sample","Length_percent_jk")],
                                     subfamily~Sample,value.var = "Length_percent_jk")
rownames(subfamily_H3K27ac_enrichment) = subfamily_H3K27ac_enrichment$subfamily
subfamily_H3K27ac_enrichment = subfamily_H3K27ac_enrichment[,2:99]
subfamily_H3K27ac_enrichment = subfamily_H3K27ac_enrichment[which(apply(subfamily_H3K27ac_enrichment,1,var) != 0),]

# PCA
sample_chromHMM_pca = prcomp(t(sample_chromHMM_3D),scale=TRUE,center=TRUE)
sample_WGBS_pca = prcomp(t(sample_WGBS_3D),scale=TRUE,center=TRUE)
subfamily_DNase_pca = prcomp(t(subfamily_DNase_enrichment),scale=TRUE,center=TRUE)
subfamily_H3K27ac_pca = prcomp(t(subfamily_H3K27ac_enrichment),scale=TRUE,center=TRUE)

# tSNE
## Samples
## chromHMM
sample_chromHMM_tsne = run_tsne(t(sample_chromHMM_3D),perplex = 10)
sample_chromHMM_tsne = merge(sample_chromHMM_tsne,metadata[,c("Sample",sample_categories)],by.x="object",by.y="Sample",all.x=TRUE)

## WGBS
sample_WGBS_tsne = run_tsne(t(sample_WGBS_3D),perplex = 10)
sample_WGBS_tsne = merge(sample_WGBS_tsne,metadata[,c("Sample",sample_categories)],by.x="object",by.y="Sample",all.x=TRUE)

## DNase
sample_DNase_tsne = run_tsne(t(subfamily_DNase_enrichment),perplex = 10)
sample_DNase_tsne = merge(sample_DNase_tsne,metadata[,c("Sample",sample_categories)],by.x="object",by.y="Sample",all.x=TRUE)
sample_DNase_tsne = merge(sample_DNase_tsne,DNase_stats[,c("Sample","Peaks")],by.x="object",by.y="Sample",all.x=TRUE)

## H3K27ac
sample_H3K27ac_tsne = run_tsne(t(subfamily_H3K27ac_enrichment),perplex = 10)
sample_H3K27ac_tsne = merge(sample_H3K27ac_tsne,metadata[,c("Sample",sample_categories)],by.x="object",by.y="Sample",all.x=TRUE)
sample_H3K27ac_tsne = merge(sample_H3K27ac_tsne,H3K27ac_stats[,c("Sample","Peaks")],by.x="object",by.y="Sample",all.x=TRUE)

## Subfamilies
## chromHMM
subfamily_chromHMM_tsne = run_tsne(t(subfamily_chromHMM_3D),perplex = 30, the = 0.5)
subfamily_chromHMM_tsne = merge(subfamily_chromHMM_tsne,rmsk_TE_subfamily[,c("subfamily","family","class_update")],by.x="object",by.y="subfamily",all.x=TRUE)

## WGBS
subfamily_WGBS_tsne = run_tsne(t(subfamily_WGBS_3D),perplex = 30, the = 0.5)
subfamily_WGBS_tsne = merge(subfamily_WGBS_tsne,rmsk_TE_subfamily[,c("subfamily","family","class_update")],by.x="object",by.y="subfamily",all.x=TRUE)

## DNase
subfamily_DNase_tsne = run_tsne(subfamily_DNase_enrichment,perplex = 30, the = 0.5)
subfamily_DNase_tsne = merge(subfamily_DNase_tsne,rmsk_TE_subfamily[,c("subfamily","family","class_update")],by.x="object",by.y="subfamily",all.x=TRUE)

## H3K27ac
subfamily_H3K27ac_tsne = run_tsne(subfamily_H3K27ac_enrichment,perplex = 30, the = 0.5)
subfamily_H3K27ac_tsne = merge(subfamily_H3K27ac_tsne,rmsk_TE_subfamily[,c("subfamily","family","class_update")],by.x="object",by.y="subfamily",all.x=TRUE)