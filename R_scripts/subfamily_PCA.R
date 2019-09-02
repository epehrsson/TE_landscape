# PCA on LOR enrichment for sample x subfamily x state combinations, by sample and subfamily

## sample_chromHMM_3D, sample_WGBS_3D, sample_DNase_enrichment, sample_H3K27ac_enrichment - matrices of LOR enrichment
## for each subfamily x state combination, by sample 

## subfamily_chromHMM_3D, subfamily_WGBS_3D, subfamily_DNase_enrichment, subfamily_H3K27ac_enrichment - matrices of LOR enrichment
## for each sample x state combination, by subfamily

## sample_chromHMM_pca, sample_WGBS_pca, sample_DNase_pca, sample_H3K27ac_pca - 
## PCA on samples, with LOR enrichment in subfamily x state combinations as variables

## subfamily_chromHMM_pca, subfamily_WGBS_pca, subfamily_DNase_pca, subfamily_H3K27ac_pca - 
## PCA on subfamilies, with LOR enrichment in sample x state combinations as variables

# Matrix of LOR enrichment of each subfamily x state combination by sample, chromHMM states
# Only subfamilies with > 30 members total and combinations with variation across samples
sample_chromHMM_3D = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% chromHMM_states & 
                                                                   rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count > THRESHOLD_IK_MEMBER),
                                                           c("subfamily","State","Sample","Enrichment")],
                                      subfamily+State~Sample,value.var = "Enrichment")
rownames(sample_chromHMM_3D) = paste(sample_chromHMM_3D$subfamily,sample_chromHMM_3D$State,sep="_")
sample_chromHMM_3D = sample_chromHMM_3D[,3:129]
sample_chromHMM_3D = sample_chromHMM_3D[which(apply(sample_chromHMM_3D,1,var) != 0),]

# Matrix of LOR enrichment of each sample x state combination by subfamily, chromHMM states
# Only subfamilies with > 30 members total and combinations with variation across subfamilies
subfamily_chromHMM_3D = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% chromHMM_states & 
                                                                      rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count > THRESHOLD_IK_MEMBER),
                                                              c("subfamily","State","Sample","Enrichment")],
                           Sample+State~subfamily,value.var = "Enrichment")
rownames(subfamily_chromHMM_3D) = paste(subfamily_chromHMM_3D$Sample,subfamily_chromHMM_3D$State,sep="_")
subfamily_chromHMM_3D = subfamily_chromHMM_3D[,3:939]
subfamily_chromHMM_3D = subfamily_chromHMM_3D[which(apply(subfamily_chromHMM_3D,1,var) != 0),]

# Matrix of LOR enrichment of each subfamily x state combination by sample, methylation states
# Only subfamilies with > 30 members with CpGs total and combinations with variation across samples
sample_WGBS_3D = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% meth_states & 
                                                               rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count_CpGs > THRESHOLD_IK_MEMBER),
                                                       c("subfamily","State","Sample","Enrichment")],
                           subfamily+State~Sample,value.var = "Enrichment")
rownames(sample_WGBS_3D) = paste(sample_WGBS_3D$subfamily,sample_WGBS_3D$State,sep="_")
sample_WGBS_3D = sample_WGBS_3D[,3:39]
sample_WGBS_3D = sample_WGBS_3D[which(apply(sample_WGBS_3D,1,var) != 0),]

# Matrix of LOR enrichment of each sample x state combination by subfamily, methylation states
# Only subfamilies with > 30 members with CpGs total and combinations with variation across subfamilies
subfamily_WGBS_3D = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State %in% meth_states & 
                                                                  rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count_CpGs > THRESHOLD_IK_MEMBER),
                                                          c("subfamily","State","Sample","Enrichment")],
                              Sample+State~subfamily,value.var = "Enrichment")
rownames(subfamily_WGBS_3D) = paste(subfamily_WGBS_3D$Sample,subfamily_WGBS_3D$State,sep="_")
subfamily_WGBS_3D = subfamily_WGBS_3D[,3:916]
subfamily_WGBS_3D = subfamily_WGBS_3D[which(apply(subfamily_WGBS_3D,1,var) != 0),]

# Matrix of LOR enrichment of each subfamily by sample, DHS
# Only subfamilies with > 30 members total and variation across samples
sample_DNase_enrichment = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State == "DNase" & 
                                                                        rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count > THRESHOLD_IK_MEMBER),
                                                                c("subfamily","Sample","Enrichment")],
                                subfamily~Sample,value.var = "Enrichment")
rownames(sample_DNase_enrichment) = sample_DNase_enrichment$subfamily
sample_DNase_enrichment = sample_DNase_enrichment[,2:54]
sample_DNase_enrichment = sample_DNase_enrichment[which(apply(sample_DNase_enrichment,1,var) != 0),]

# Matrix of LOR enrichment of each sample by subfamily, DHS
# Only subfamilies with > 30 members total and samples with variation across subfamilies
subfamily_DNase_enrichment = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State == "DNase" & 
                                                                           rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count > THRESHOLD_IK_MEMBER),
                                                                   c("subfamily","Sample","Enrichment")],
                                   Sample~subfamily,value.var = "Enrichment")
rownames(subfamily_DNase_enrichment) = subfamily_DNase_enrichment$Sample
subfamily_DNase_enrichment = subfamily_DNase_enrichment[,2:938]
subfamily_DNase_enrichment = subfamily_DNase_enrichment[which(apply(subfamily_DNase_enrichment,1,var) != 0),]

# Matrix of LOR enrichment of each subfamily by sample, H3K27ac
# Only subfamilies with > 30 members total and variation across samples
sample_H3K27ac_enrichment = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State == "H3K27ac" & 
                                                                          rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count > THRESHOLD_IK_MEMBER),
                                                                  c("subfamily","Sample","Enrichment")],
                                subfamily~Sample,value.var = "Enrichment")
rownames(sample_H3K27ac_enrichment) = sample_H3K27ac_enrichment$subfamily
sample_H3K27ac_enrichment = sample_H3K27ac_enrichment[,2:99]
sample_H3K27ac_enrichment = sample_H3K27ac_enrichment[which(apply(sample_H3K27ac_enrichment,1,var) != 0),]

# Matrix of LOR enrichment of each sample by subfamily, H3K27ac
# Only subfamilies with > 30 members total and samples with variation across subfamilies
subfamily_H3K27ac_enrichment = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State == "H3K27ac" & 
                                                                             rmsk_TE_subfamily[match(subfamily_state_sample_combined$subfamily,rmsk_TE_subfamily$subfamily),]$Count > THRESHOLD_IK_MEMBER),
                                                                     c("subfamily","Sample","Enrichment")],
                                     Sample~subfamily,value.var = "Enrichment")
rownames(subfamily_H3K27ac_enrichment) = subfamily_H3K27ac_enrichment$Sample
subfamily_H3K27ac_enrichment = subfamily_H3K27ac_enrichment[,2:938]
subfamily_H3K27ac_enrichment = subfamily_H3K27ac_enrichment[which(apply(subfamily_H3K27ac_enrichment,1,var) != 0),]

# PCA on samples, with LOR enrichment in subfamily x state combinations as variables
# scaled and centered, with eigenvalues/eigenvectors calculated

## chromHMM
sample_chromHMM_pca = prcomp(t(sample_chromHMM_3D),scale=TRUE,center=TRUE)
sample_chromHMM_pca = format_pca(sample_chromHMM_pca,metadata,"Sample")

## Methylation
sample_WGBS_pca = prcomp(t(sample_WGBS_3D),scale=TRUE,center=TRUE)
sample_WGBS_pca = format_pca(sample_WGBS_pca,metadata,"Sample")

## DHS
sample_DNase_pca = prcomp(t(sample_DNase_enrichment),scale=TRUE,center=TRUE)
sample_DNase_pca = format_pca(sample_DNase_pca,metadata,"Sample")

## H3K27ac
sample_H3K27ac_pca = prcomp(t(sample_H3K27ac_enrichment),scale=TRUE,center=TRUE)
sample_H3K27ac_pca = format_pca(sample_H3K27ac_pca,metadata,"Sample")

# PCA on subfamilies, with LOR enrichment in sample x state combinations as variables
# scaled and centered, with eigenvalues/eigenvectors calculated

## chromHMM
subfamily_chromHMM_pca = prcomp(t(subfamily_chromHMM_3D),scale=TRUE,center=TRUE)
subfamily_chromHMM_pca = format_pca(subfamily_chromHMM_pca,rmsk_TE_subfamily,"subfamily")

## Methylation
subfamily_WGBS_pca = prcomp(t(subfamily_WGBS_3D),scale=TRUE,center=TRUE)
subfamily_WGBS_pca = format_pca(subfamily_WGBS_pca,rmsk_TE_subfamily,"subfamily")

## DHS
subfamily_DNase_pca = prcomp(t(subfamily_DNase_enrichment),scale=TRUE,center=TRUE)
subfamily_DNase_pca = format_pca(subfamily_DNase_pca,rmsk_TE_subfamily,"subfamily")

## H3K27ac
subfamily_H3K27ac_pca = prcomp(t(subfamily_H3K27ac_enrichment),scale=TRUE,center=TRUE)
subfamily_H3K27ac_pca = format_pca(subfamily_H3K27ac_pca,rmsk_TE_subfamily,"subfamily")