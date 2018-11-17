# Perform PCA at the subfamily level

# chromHMM 3D clustering - by sample
sample_chromHMM_3D = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% chromHMM_states & 
                                                                   rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count > THRESHOLD_IK_MEMBER),
                                                           c("subfamily","State","Sample","Enrichment")],
                                      subfamily+State~Sample,value.var = "Enrichment")
rownames(sample_chromHMM_3D) = paste(sample_chromHMM_3D$subfamily,sample_chromHMM_3D$State,sep="_")
sample_chromHMM_3D = sample_chromHMM_3D[,3:129]
sample_chromHMM_3D = sample_chromHMM_3D[which(apply(sample_chromHMM_3D,1,var) != 0),]

# chromHMM 3D clustering - by subfamily
subfamily_chromHMM_3D = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% chromHMM_states & 
                                                                      rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count > THRESHOLD_IK_MEMBER),
                                                              c("subfamily","State","Sample","Enrichment")],
                           Sample+State~subfamily,value.var = "Enrichment")
rownames(subfamily_chromHMM_3D) = paste(subfamily_chromHMM_3D$Sample,subfamily_chromHMM_3D$State,sep="_")
subfamily_chromHMM_3D = subfamily_chromHMM_3D[,3:939]
subfamily_chromHMM_3D = subfamily_chromHMM_3D[which(apply(subfamily_chromHMM_3D,1,var) != 0),]

# WGBS 3D clustering - by sample
sample_WGBS_3D = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% meth_states & 
                                                               rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count_CpGs > THRESHOLD_IK_MEMBER),
                                                       c("subfamily","State","Sample","Enrichment")],
                           subfamily+State~Sample,value.var = "Enrichment")
rownames(sample_WGBS_3D) = paste(sample_WGBS_3D$subfamily,sample_WGBS_3D$State,sep="_")
sample_WGBS_3D = sample_WGBS_3D[,3:39]
sample_WGBS_3D = sample_WGBS_3D[which(apply(sample_WGBS_3D,1,var) != 0),]

# WGBS 3D clustering - by subfamily
subfamily_WGBS_3D = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% meth_states & 
                                                                  rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count_CpGs > THRESHOLD_IK_MEMBER),
                                                          c("subfamily","State","Sample","Enrichment")],
                              Sample+State~subfamily,value.var = "Enrichment")
rownames(subfamily_WGBS_3D) = paste(subfamily_WGBS_3D$Sample,subfamily_WGBS_3D$State,sep="_")
subfamily_WGBS_3D = subfamily_WGBS_3D[,3:916]
subfamily_WGBS_3D = subfamily_WGBS_3D[which(apply(subfamily_WGBS_3D,1,var) != 0),]

## DNase - by sample
sample_DNase_enrichment = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State == "DNase" & 
                                                                        rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count > THRESHOLD_IK_MEMBER),
                                                                c("subfamily","Sample","Enrichment")],
                                subfamily~Sample,value.var = "Enrichment")
rownames(sample_DNase_enrichment) = sample_DNase_enrichment$subfamily
sample_DNase_enrichment = sample_DNase_enrichment[,2:54]
sample_DNase_enrichment = sample_DNase_enrichment[which(apply(sample_DNase_enrichment,1,var) != 0),]

## DNase - by subfamily
subfamily_DNase_enrichment = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State == "DNase" & 
                                                                           rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count > THRESHOLD_IK_MEMBER),
                                                                   c("subfamily","Sample","Enrichment")],
                                   Sample~subfamily,value.var = "Enrichment")
rownames(subfamily_DNase_enrichment) = subfamily_DNase_enrichment$Sample
subfamily_DNase_enrichment = subfamily_DNase_enrichment[,2:938]
subfamily_DNase_enrichment = subfamily_DNase_enrichment[which(apply(subfamily_DNase_enrichment,1,var) != 0),]

## H3K27ac - by sample
sample_H3K27ac_enrichment = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State == "H3K27ac" & 
                                                                          rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count > THRESHOLD_IK_MEMBER),
                                                                  c("subfamily","Sample","Enrichment")],
                                subfamily~Sample,value.var = "Enrichment")
rownames(sample_H3K27ac_enrichment) = sample_H3K27ac_enrichment$subfamily
sample_H3K27ac_enrichment = sample_H3K27ac_enrichment[,2:99]
sample_H3K27ac_enrichment = sample_H3K27ac_enrichment[which(apply(sample_H3K27ac_enrichment,1,var) != 0),]

## H3K27ac - by subfamily
subfamily_H3K27ac_enrichment = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State == "H3K27ac" & 
                                                                             rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count > THRESHOLD_IK_MEMBER),
                                                                     c("subfamily","Sample","Enrichment")],
                                     Sample~subfamily,value.var = "Enrichment")
rownames(subfamily_H3K27ac_enrichment) = subfamily_H3K27ac_enrichment$Sample
subfamily_H3K27ac_enrichment = subfamily_H3K27ac_enrichment[,2:938]
subfamily_H3K27ac_enrichment = subfamily_H3K27ac_enrichment[which(apply(subfamily_H3K27ac_enrichment,1,var) != 0),]

# PCA
## Samples
sample_chromHMM_pca = prcomp(t(sample_chromHMM_3D),scale=TRUE,center=TRUE)
sample_chromHMM_pca = format_pca(sample_chromHMM_pca,metadata,"Sample")

sample_WGBS_pca = prcomp(t(sample_WGBS_3D),scale=TRUE,center=TRUE)
sample_WGBS_pca = format_pca(sample_WGBS_pca,metadata,"Sample")

sample_DNase_pca = prcomp(t(sample_DNase_enrichment),scale=TRUE,center=TRUE)
sample_DNase_pca = format_pca(sample_DNase_pca,metadata,"Sample")

sample_H3K27ac_pca = prcomp(t(sample_H3K27ac_enrichment),scale=TRUE,center=TRUE)
sample_H3K27ac_pca = format_pca(sample_H3K27ac_pca,metadata,"Sample")

## Subfamily
subfamily_chromHMM_pca = prcomp(t(subfamily_chromHMM_3D),scale=TRUE,center=TRUE)
subfamily_chromHMM_pca = format_pca(subfamily_chromHMM_pca,rmsk_TE_subfamily,"subfamily")

subfamily_WGBS_pca = prcomp(t(subfamily_WGBS_3D),scale=TRUE,center=TRUE)
subfamily_WGBS_pca = format_pca(subfamily_WGBS_pca,rmsk_TE_subfamily,"subfamily")

subfamily_DNase_pca = prcomp(t(subfamily_DNase_enrichment),scale=TRUE,center=TRUE)
subfamily_DNase_pca = format_pca(subfamily_DNase_pca,rmsk_TE_subfamily,"subfamily")

subfamily_H3K27ac_pca = prcomp(t(subfamily_H3K27ac_enrichment),scale=TRUE,center=TRUE)
subfamily_H3K27ac_pca = format_pca(subfamily_H3K27ac_pca,rmsk_TE_subfamily,"subfamily")
