# Read in human-mouse orthologous TEs and metadata

# Corresponding human-mouse samples
human_mouse_samples = read.table("sample_lists/human_mouse_samples.txt",sep='\t')
colnames(human_mouse_samples) = c("Human_sample","Mouse_sample_WGBS","Mouse_sample_chromHMM")
human_mouse_samples = merge(human_mouse_samples,metadata[,c("Sample",sample_categories)],by.x="Human_sample",by.y="Sample")

# Expanded human-mouse samples (chromHMM)
a = expand.grid(Sample1=human_mouse_samples$Human_sample,Sample2=human_mouse_samples$Mouse_sample_chromHMM)
a = merge(a,human_mouse_samples[,c("Human_sample","Mouse_sample_chromHMM")],by.x=c("Sample1","Sample2"),by.y=c("Human_sample","Mouse_sample_chromHMM"),all=TRUE)
a$Test = rep("Paired",dim(a)[1])
a = merge(a,human_mouse_samples[,c("Human_sample",sample_categories)],by.x="Sample1",by.y="Human_sample")
a = merge(a,human_mouse_samples[,c("Mouse_sample_chromHMM",sample_categories)],by.x="Sample2",by.y="Mouse_sample_chromHMM")

b = expand.grid(Sample1=human_mouse_samples$Mouse_sample_chromHMM,Sample2=human_mouse_samples$Mouse_sample_chromHMM)
b$Test = rep("Mouse",dim(b)[1])
b = merge(b,human_mouse_samples[,c("Mouse_sample_chromHMM",sample_categories)],by.x="Sample1",by.y="Mouse_sample_chromHMM")
b = merge(b,human_mouse_samples[,c("Mouse_sample_chromHMM",sample_categories)],by.x="Sample2",by.y="Mouse_sample_chromHMM")

c = expand.grid(Sample1=human_mouse_samples$Human_sample,Sample2=human_mouse_samples$Human_sample)
c$Test = rep("Human",dim(c)[1])
c = merge(c,human_mouse_samples[,c("Human_sample",sample_categories)],by.x="Sample1",by.y="Human_sample")
c = merge(c,human_mouse_samples[,c("Human_sample",sample_categories)],by.x="Sample2",by.y="Human_sample")

human_mouse_samples_expand = rbind(a,b,c)
rm(list=c("a","b","c"))
human_mouse_samples_expand = human_mouse_samples_expand[which(human_mouse_samples_expand$Sample1 != human_mouse_samples_expand$Sample2),]
human_mouse_samples_expand$Tissue = ifelse(human_mouse_samples_expand$Anatomy.x == human_mouse_samples_expand$Anatomy.y,"Same",
                                           ifelse(human_mouse_samples_expand$Anatomy.x %in% c("GI_INTESTINE","GI_COLON") & 
                                                    human_mouse_samples_expand$Anatomy.y %in% c("GI_INTESTINE","GI_COLON"),"Same","Different"))
human_mouse_samples_expand$Age = ifelse(human_mouse_samples_expand$Age.x == human_mouse_samples_expand$Age.y,"Same","Different")
human_mouse_samples_expand$Tissue.Age = factor(apply(human_mouse_samples_expand,1,function(x) paste(x[16],x[17],sep="_")),levels=c("Same_Same","Same_Different","Different_Same","Different_Different"))

# Expanded human-mouse samples (WGBS)
human_mouse_samples_WGBS = human_mouse_samples_expand[which(human_mouse_samples_expand$Sample1 %in% c(as.vector(na.omit(human_mouse_samples)$Human_sample),as.vector(na.omit(human_mouse_samples)$Mouse_sample_chromHMM)) & 
                                                              human_mouse_samples_expand$Sample2 %in% c(as.vector(na.omit(human_mouse_samples)$Human_sample),as.vector(na.omit(human_mouse_samples)$Mouse_sample_chromHMM))),] 
human_mouse_samples_WGBS$Sample1 = mapvalues(human_mouse_samples_WGBS$Sample1,as.vector(na.omit(human_mouse_samples)$Mouse_sample_chromHMM),as.vector(na.omit(human_mouse_samples)$Mouse_sample_WGBS))
human_mouse_samples_WGBS$Sample2 = mapvalues(human_mouse_samples_WGBS$Sample2,as.vector(na.omit(human_mouse_samples)$Mouse_sample_chromHMM),as.vector(na.omit(human_mouse_samples)$Mouse_sample_WGBS))

# Human TEs with mouse orthologs (mm10)
human_mouse_orthologs_mm10 = read.table("Mouse/liftover/hg19_mm10_TE_intersect_same.bed",sep='\t')
colnames(human_mouse_orthologs_mm10) = c("human_chr_mm10","human_start_mm10","human_stop_mm10","human_strand_mm10",hg19_coordinates,mm10_coordinates,"overlap")
human_mouse_orthologs_mm10$class_update = convert_class(human_mouse_orthologs_mm10$human_class)
human_mouse_orthologs_mm10$human_TE_mm10 = apply(human_mouse_orthologs_mm10,1,function(x) paste(x[1],x[2],x[3],x[4],sep="_"))
human_mouse_orthologs_mm10$human_TE_hg19 = apply(human_mouse_orthologs_mm10,1,function(x) paste(x[5],x[6],x[7],x[8],x[9],x[10],x[11],sep="_"))
human_mouse_orthologs_mm10$mouse_TE_mm10 = apply(human_mouse_orthologs_mm10,1,function(x) paste(x[12],x[13],x[14],x[15],x[16],x[17],x[18],sep="_"))
human_mouse_orthologs_mm10 = merge(human_mouse_orthologs_mm10,rmsk_TE[,c(TE_coordinates,"JC_distance")],by.x=hg19_coordinates,by.y=TE_coordinates[c(1:4,6,5,7)])

# Mouse TEs (mm10)
mm10_rmsk_TE = read.table("features/mouse/mm10_rmsk_TE.txt",sep='\t')
colnames(mm10_rmsk_TE) = TE_coordinates[c(1:4,6,5,7)]