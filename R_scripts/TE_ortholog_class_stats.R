# Add human-mouse ortholog stats to TE class stats

# Mouse orthologs and subfamilies by class
rmsk_TE_class$Mouse_orthologs = table(human_mouse_orthologs_mm10$class_update)[as.vector(rmsk_TE_class$class_update)]
rmsk_TE_class = merge(rmsk_TE_class,ddply(rmsk_TE_subfamily[which(rmsk_TE_subfamily$subfamily %in% mm10_rmsk_TE$subfamily),],~class_update,summarize,
                                          Mouse_ortholog_subfamily = length(subfamily)),by="class_update",all.x=TRUE)
rmsk_TE_class[which(is.na(rmsk_TE_class$Mouse_ortholog_subfamily)),]$Mouse_ortholog_subfamily = 0