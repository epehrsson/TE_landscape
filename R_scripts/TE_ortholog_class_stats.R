# Add human-mouse ortholog stats to TE class stats

#source("R_scripts/TE_class_stats.R")
#source("R_scripts/TE_subfamily_stats.R")
#source("R_scripts/human_mouse_ortholog_WGBS.R")

# Mouse orthologs and subfamilies by class
test = table(unique(human_mouse_orthologs_mm10[,5:11])$human_class)
test[12] = sum(test[c(2,4,6:7,9:11)])
names(test)[12] = "Other"
rmsk_TE_class$Mouse_orthologs = test[as.vector(rmsk_TE_class$class_update)]
rmsk_TE_class[which(is.na(rmsk_TE_class$Mouse_orthologs)),]$Mouse_orthologs = 0
rm(test)

mm10_rmsk_TE = read.table("features/mouse/TEs/mm10_rmsk_TE.txt",sep='\t')
colnames(mm10_rmsk_TE) = c("chromosome","start","stop","subfamily","class","family","strand")

rmsk_TE_class = merge(rmsk_TE_class,ddply(rmsk_TE_subfamily[which(rmsk_TE_subfamily$subfamily %in% mm10_rmsk_TE$subfamily),],~class_update,summarize,Mouse_ortholog_subfamily = length(subfamily)),by="class_update",all.x=TRUE)
rmsk_TE_class[which(is.na(rmsk_TE_class$Mouse_ortholog_subfamily)),]$Mouse_ortholog_subfamily = 0
