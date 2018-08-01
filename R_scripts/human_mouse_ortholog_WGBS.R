# Comparison of human-mouse TE ortholog methylation state
# See 5/11/2017, 5/12/2017, 5/22/2017, 5/29/2017, 5/30/2017, 7/21/2017, 7/26/2017, 7/31/2017, 8/1/2017

#load("R_datasets/TE_meth_average.RData")
#source("R_scripts/WGBS_subfamily_enrichment_TE.R")

# Corresponding human-mouse samples
human_mouse_samples = read.table("Mouse/human_mouse_samples.txt",sep='\t',header=TRUE)
colnames(human_mouse_samples) = c("Human_sample","Mouse_sample_WGBS","Mouse_sample_chromHMM")
human_mouse_samples$Test = rep("Paired",dim(human_mouse_samples)[1])
human_mouse_samples = merge(human_mouse_samples,metadata[,c("Sample",sample_categories)],by.x="Human",by.y="Sample")

# Expanded human-mouse samples
a = expand.grid(unique(human_mouse_samples$Mouse),unique(human_mouse_samples$Human))
a = merge(a,human_mouse_samples[,c("Mouse","Human","Test")],by.x=c("Var1","Var2"),by.y=c("Mouse","Human"),all=TRUE)
a[is.na(a)] = "Random_Pair"
b = expand.grid(unique(human_mouse_samples$Mouse),unique(human_mouse_samples$Mouse))
b$Test = rep("Mouse",dim(b)[1])
c = expand.grid(unique(human_mouse_samples$Human),unique(human_mouse_samples$Human))
c$Test = rep("Human",dim(c)[1])
human_mouse_samples_expand = rbind(a,b,c)
rm(list=c("a","b","c"))
colnames(human_mouse_samples_expand)[1:2] = c("Sample1","Sample2")
human_mouse_samples_expand = human_mouse_samples_expand[which(human_mouse_samples_expand$Sample1 != human_mouse_samples_expand$Sample2),]
human_mouse_samples_expand = merge(human_mouse_samples_expand,unique(melt(human_mouse_samples,id.vars=c("Test","Anatomy")))[1:22,c("value","Anatomy")],by.x="Sample1",by.y="value")
human_mouse_samples_expand = merge(human_mouse_samples_expand,unique(melt(human_mouse_samples,id.vars=c("Test","Anatomy")))[1:22,c("value","Anatomy")],by.x="Sample2",by.y="value")
human_mouse_samples_expand$Tissue = ifelse(human_mouse_samples_expand$Anatomy.x == human_mouse_samples_expand$Anatomy.y,"Same",
                                           ifelse(human_mouse_samples_expand$Anatomy.x %in% c("GI_STOMACH","GI_INTESTINE","GI_COLON") & 
                                                    human_mouse_samples_expand$Anatomy.y %in% c("GI_STOMACH","GI_INTESTINE","GI_COLON"),"Same","Different"))
  
# Human TEs with mouse orthologs (mm10)
human_mouse_orthologs_mm10 = read.table("Mouse/liftover/hg19_mm10_TE_intersect_same.bed",sep='\t')
colnames(human_mouse_orthologs_mm10) = c("human_chr_mm10","human_start_mm10","human_stop_mm10","human_strand_mm10",hg19_coordinates,mm10_coordinates,"overlap")
human_mouse_orthologs_mm10$class_update = convert_class(human_mouse_orthologs_mm10$human_class)
human_mouse_orthologs_mm10$human_TE_mm10 = apply(human_mouse_orthologs_mm10,1,function(x) paste(x[1],x[2],x[3],x[4],sep="_"))
human_mouse_orthologs_mm10$human_TE_hg19 = apply(human_mouse_orthologs_mm10,1,function(x) paste(x[5],x[6],x[7],x[8],x[9],x[10],x[11],sep="_"))
human_mouse_orthologs_mm10$mouse_TE_mm10 = apply(human_mouse_orthologs_mm10,1,function(x) paste(x[12],x[13],x[14],x[15],x[16],x[17],x[18],sep="_"))

# Mouse TEs (mm10)
mm10_rmsk_TE = read.table("features/mouse/mm10_rmsk_TE.txt",sep='\t')
colnames(mm10_rmsk_TE) = TE_coordinates[c(1:4,6,5,7)]

# mm10 TEs with WGBS average methylation
mm10_rmsk_TE_WGBS = read.table("Mouse/WGBS/mm10_rmsk_TE_WGBS_avg.txt",sep='\t',header=TRUE)
mm10_rmsk_TE_WGBS[,8:16] = mm10_rmsk_TE_WGBS[,8:16]/100
mm10_rmsk_TE_CpG_count = read.table("Mouse/WGBS/mm10_rmsk_TE_CpG_count.txt",sep='\t')
colnames(mm10_rmsk_TE_CpG_count) = c(TE_coordinates[c(1:4,6,5,7)],"CpGs")
mm10_rmsk_TE_CpG_count$CpGs = mm10_rmsk_TE_CpG_count$CpGs/2
mm10_rmsk_TE_WGBS = merge(mm10_rmsk_TE_WGBS,mm10_rmsk_TE_CpG_count,by=TE_coordinates)
rm(mm10_rmsk_TE_CpG_count)

# Orthologous TEs with methylation level (only those with CpGs in both)
human_mouse_orthologs_WGBS = merge(human_mouse_orthologs_mm10[,c("human_chr_mm10","human_start_mm10","human_stop_mm10","human_strand_mm10",hg19_coordinates,mm10_coordinates,"overlap")],
                                   TE_meth_average[,c(TE_coordinates,as.vector(human_mouse_samples$Human))],by.y=TE_coordinates[c(1:4,6,5,7)],by.x=hg19_coordinates)
human_mouse_orthologs_WGBS = merge(human_mouse_orthologs_WGBS,mm10_rmsk_TE_WGBS,by.y=TE_coordinates[c(1:4,6,5,7)],by.x=mm10_coordinates)

# Paired samples only
human_mouse_orthologs_WGBS_paired = melt(human_mouse_orthologs_WGBS,id.vars=c("human_chr_mm10","human_start_mm10","human_stop_mm10","human_strand_mm10",hg19_coordinates,mm10_coordinates,"overlap","CpGs",
                                                                              as.vector(human_mouse_samples$Human)))
colnames(human_mouse_orthologs_WGBS_paired)[34:35] = c("Mouse_sample","Mouse_methylation")
human_mouse_orthologs_WGBS_paired = melt(human_mouse_orthologs_WGBS_paired,id.vars=c("human_chr_mm10","human_start_mm10","human_stop_mm10","human_strand_mm10",hg19_coordinates,mm10_coordinates,"overlap","CpGs","Mouse_sample","Mouse_methylation"))
colnames(human_mouse_orthologs_WGBS_paired)[23:24] = c("Human_sample","Human_methylation")
human_mouse_orthologs_WGBS_paired = merge(human_mouse_orthologs_WGBS_paired,human_mouse_samples[,c("Mouse","Human")],by.x=c("Mouse_sample","Human_sample"),by.y=c("Mouse","Human"))

# Adding methylation state
human_mouse_orthologs_WGBS_paired$Human_state = factor(ifelse(is.na(human_mouse_orthologs_WGBS_paired$Human_methylation),"Missing",
                                                       ifelse(human_mouse_orthologs_WGBS_paired$Human_methylation > 0.7,"Hypermethylated",
                                                       ifelse(human_mouse_orthologs_WGBS_paired$Human_methylation < 0.3,"Hypomethylated",
                                                              ifelse(human_mouse_orthologs_WGBS_paired$Human_methylation <= 0.7 & human_mouse_orthologs_WGBS_paired$Human_methylation >= 0.3,"Intermediate",NA)))),
                                                       levels=meth_states)
human_mouse_orthologs_WGBS_paired$Mouse_state = factor(ifelse(is.na(human_mouse_orthologs_WGBS_paired$Mouse_methylation),"Missing",
                                                       ifelse(human_mouse_orthologs_WGBS_paired$Mouse_methylation > 0.7,"Hypermethylated",
                                                       ifelse(human_mouse_orthologs_WGBS_paired$Mouse_methylation < 0.3,"Hypomethylated",
                                                              ifelse(human_mouse_orthologs_WGBS_paired$Mouse_methylation <= 0.7 & human_mouse_orthologs_WGBS_paired$Mouse_methylation >= 0.3,"Intermediate",NA)))),
                                                       levels=meth_states)

# Adding chromHMM sample
human_mouse_orthologs_WGBS_paired$Mouse_sample_chromHMM = ifelse(human_mouse_orthologs_WGBS_paired$Mouse_sample %in% c("ENCFF592SYK","ENCFF796POX","ENCFF446IJS"),"ENCFF640YBG",
                                                          ifelse(human_mouse_orthologs_WGBS_paired$Mouse_sample %in% c("ENCFF895RRK","ENCFF644XMA"),"ENCFF715UGP",
                                                                 ifelse(human_mouse_orthologs_WGBS_paired$Mouse_sample == "ENCFF467UEZ","ENCFF880AMI",
                                                                        ifelse(human_mouse_orthologs_WGBS_paired$Mouse_sample == "ENCFF770AUO","ENCFF539POM",NA))))
# Adding age
human_mouse_orthologs_WGBS_paired = merge(human_mouse_orthologs_WGBS_paired,rmsk_TE[,c(TE_coordinates,"JC_distance")],by.x=hg19_coordinates,by.y=TE_coordinates[c(1:4,6,5,7)])

# TE orthologs hypomethylated in both human and mouse
human_mouse_orthologs_WGBS_hypo = human_mouse_orthologs_WGBS_paired[which(human_mouse_orthologs_WGBS_paired$Mouse_state == "Hypomethylated" & human_mouse_orthologs_WGBS_paired$Human_state == "Hypomethylated"),]

# Spearman correlation between all TEs for each sample pair
human_mouse_spearman = apply(human_mouse_orthologs_WGBS[,c(unique(as.vector(human_mouse_samples$Mouse)),as.vector(human_mouse_samples$Human))],2,function(x) 
  apply(human_mouse_orthologs_WGBS[,c(unique(as.vector(human_mouse_samples$Mouse)),as.vector(human_mouse_samples$Human))],2,function(y) unlist(cor.test(x,y,method="spearman"))["estimate.rho"]))
human_mouse_spearman = melt(as.matrix(human_mouse_spearman))
colnames(human_mouse_spearman) = c("Sample1","Sample2","Correlation")
human_mouse_spearman = merge(human_mouse_spearman,human_mouse_samples_expand)
human_mouse_spearman$Correlation = as.numeric(as.character(human_mouse_spearman$Correlation))