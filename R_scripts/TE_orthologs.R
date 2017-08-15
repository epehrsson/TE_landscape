# Human-mouse ortholog coordinates
# See 9/28/2016, 2/9/2017

# Human-mouse orthologous TEs
human_mouse_orthologs = read.table("Mouse/hg19_mm9_TE_intersect_same.bed",sep='\t')
colnames(human_mouse_orthologs) = c("mm9_chr_human","mm9_start_human","mm9_stop_human","mm9_strand_human","hg19_chr_human","hg19_start_human","hg19_stop_human","Subfamily_human","Class_human","Family_human","hg19_strand_human","mm9_chr_mouse","mm9_start_mouse","mm9_stop_mouse","Subfamily_mouse","Class_mouse","Family_mouse","mm9_strand_mouse","Overlap")

# Human TEs with mouse orthologs
human_mouse_orthologs_other = read.table("Mouse/hg19_mm9_other_intersect_same.bed",sep='\t')
human_mouse_orthologs_TE_other = read.table("Mouse/hg19_mm9_TE_other_intersect_same.bed",sep='\t')
colnames(human_mouse_orthologs_other) = c("mm9_chr_human","mm9_start_human","mm9_stop_human","mm9_strand_human","hg19_chr_human","hg19_start_human","hg19_stop_human","Subfamily_human","Class_human","Family_human","hg19_strand_human","mm9_chr_mouse","mm9_start_mouse","mm9_stop_mouse","Subfamily_mouse","Class_mouse","Family_mouse","mm9_strand_mouse","Overlap")
colnames(human_mouse_orthologs_TE_other) = c("mm9_chr_human","mm9_start_human","mm9_stop_human","mm9_strand_human","hg19_chr_human","hg19_start_human","hg19_stop_human","Subfamily_human","Class_human","Family_human","hg19_strand_human","mm9_chr_mouse","mm9_start_mouse","mm9_stop_mouse","Subfamily_mouse","Class_mouse","Family_mouse","mm9_strand_mouse","Overlap")
human_mouse_orthologs_TEother = rbind(human_mouse_orthologs,human_mouse_orthologs_other,human_mouse_orthologs_TE_other)
