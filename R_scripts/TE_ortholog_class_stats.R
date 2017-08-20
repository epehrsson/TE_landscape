# Add human-mouse ortholog stats to TE class stats

source("R_scripts/TE_class_stats.R")

# Mouse orthologs and subfamilies by class
human_mouse_orthologs_mm10 = read.table("Mouse/liftover/hg19_mm10_TE_intersect_same.bed",sep='\t')
colnames(human_mouse_orthologs_mm10) = c("human_chr_mm10","human_start_mm10","human_stop_mm10","human_strand_mm10","human_chr_hg19","human_start_hg19","human_stop_hg19","human_subfamily","human_class","human_family","human_strand_hg19","mouse_chr_mm10","mouse_start_mm10","mouse_stop_mm10","mouse_subfamily","mouse_class","mouse_family","mouse_strand_mm10","overlap")

test = table(unique(human_mouse_orthologs_mm10[,5:11])$human_class)
test[12] = sum(test[c(2,4,6:7,9:11)])
names(test)[12] = "Other"
rmsk_TE_class$Mouse_orthologs = test[as.vector(rmsk_TE_class$class_update)]
rmsk_TE_class[which(is.na(rmsk_TE_class$Mouse_orthologs)),]$Mouse_orthologs = 0
rm(test)

mm10_rmsk_TE = read.table("features/mouse/TEs/mm10_rmsk_TE.txt",sep='\t')
colnames(mm10_rmsk_TE) = c("chromosome","start","stop","subfamily","class","family","strand")

source("R_scripts/TE_subfamily_stats.R")

rmsk_TE_class = merge(rmsk_TE_class,ddply(rmsk_TE_subfamily[which(rmsk_TE_subfamily$subfamily %in% mm10_rmsk_TE$subfamily),],~class_update,summarize,Mouse_ortholog_subfamily = length(subfamily)),by="class_update",all.x=TRUE)
rmsk_TE_class[which(is.na(rmsk_TE_class$Mouse_ortholog_subfamily)),]$Mouse_ortholog_subfamily = 0