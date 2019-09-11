# Creates dataframes of human and corresponding mouse samples, mm10 TEs, and hg19/mm10 orthologous TEs

## human_mouse_samples - Table of mouse ENCODE samples corresponding to human Roadmap samples
## human_mouse_orthologs_mm10 - Dataframe of human TEs and mouse (mm10) orthologs, including human TE age
## mm10_rmsk_TE - Dataframe of mouse TEs (mm10)


# Table of mouse ENCODE samples (chromHMM and WGBS) corresponding to human Roadmap samples
human_mouse_samples = read.table("sample_lists/human_mouse_samples.txt",sep='\t',col.names=c("Human_sample","Mouse_sample_WGBS","Mouse_sample_chromHMM"))
## Add sample categories for human samples
human_mouse_samples = merge(human_mouse_samples,metadata[,c("Sample",sample_categories)],by.x="Human_sample",by.y="Sample")


# Dataframe of human TEs and mouse (mm10) orthologs
human_mouse_orthologs_mm10 = read.table("Mouse/liftover/hg19_mm10_TE_intersect_same.bed",sep='\t',
                                        col.names=c("human_chr_mm10","human_start_mm10","human_stop_mm10","human_strand_mm10",hg19_coordinates,mm10_coordinates,"overlap"))
human_mouse_orthologs_mm10$class_update = convert_class(human_mouse_orthologs_mm10$human_class)
## Create columns with TE coordinates for hg19 TEs, mm10 TEs, and hg19 TEs lifted over to mm10
human_mouse_orthologs_mm10$human_TE_mm10 = apply(human_mouse_orthologs_mm10,1,function(x) paste(x[1],x[2],x[3],x[4],sep="_"))
human_mouse_orthologs_mm10$human_TE_hg19 = apply(human_mouse_orthologs_mm10,1,function(x) paste(x[5],x[6],x[7],x[8],x[9],x[10],x[11],sep="_"))
human_mouse_orthologs_mm10$mouse_TE_mm10 = apply(human_mouse_orthologs_mm10,1,function(x) paste(x[12],x[13],x[14],x[15],x[16],x[17],x[18],sep="_"))
## Add JC distance for each human TE
human_mouse_orthologs_mm10 = merge(human_mouse_orthologs_mm10,rmsk_TE[,c(TE_coordinates,"JC_distance")],by.x=hg19_coordinates,by.y=TE_coordinates[c(1:4,6,5,7)])


# Dataframe of mouse TEs (mm10)
mm10_rmsk_TE = read.table("features/mouse/mm10_rmsk_TE.txt",sep='\t',col.names=TE_coordinates[c(1:4,6,5,7)])