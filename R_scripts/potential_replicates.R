# Get metadata table
metadata_table = ldply(apply(metadata[,sample_categories],2,as.data.frame(table)))
colnames(metadata_table) = c("Category","Grouping","Samples")

# Get samples in state
active_reg_samples = melt(chromHMM_TE_state[,c(TE_coordinates,states[c(1:3,6:7)])],id.vars=TE_coordinates)
colnames(active_reg_samples)[8:9] = c("State","Total")
active_reg_samples = active_reg_samples[which(active_reg_samples$Total > 0),]

# Load TEs in state
active_reg = read.table("chromHMM/chromHMM_active_reg.txt",sep='\t')[,c(1:8,10)]
colnames(active_reg) = c(TE_coordinates[c(1:4,6,5,7)],"Sample","State")

# Add metadata
active_reg = merge(active_reg,metadata[,c("Sample",sample_categories)],by="Sample")
active_reg = melt(active_reg,id.vars=c(TE_coordinates[c(1:4,6,5,7)],"Sample","State"))
colnames(active_reg)[10:11] = c("Category","Grouping")

# Count samples in Group
active_reg_group = ddply(active_reg,.(chromosome,start,stop,subfamily,family,class,strand,State,Category,Grouping),
                                summarise,Count=length(Grouping))
active_reg_group = merge(active_reg_group,active_reg_samples,by=c(TE_coordinates,"State"))
colnames(active_reg_group)[11] = "Total"

# Add Group totals
active_reg_group = merge(active_reg_group,metadata_table,by=c("Category","Grouping"))
colnames(active_reg_group)[12] = "Category.Total"

# Filtering to TEs not in an active regulatory state outside the category
# and removing categories with only one sample
active_reg_excl = active_reg_group[which(active_reg_group$Count == active_reg_group$Total & 
                                       active_reg_group$Category.Total != 1),]

# Fraction of TEs in either pair found in both
active_reg_table = ddply(active_reg_excl,.(Category,Grouping,Category.Total),summarise,
      Single=sum(Count == 1),Both=sum(Count >= 2))
active_reg_table$Fraction = active_reg_table$Both/(active_reg_table$Single+active_reg_table$Both)


# For TEs in 2+ samples, how often is it the same category?
active_reg_multi = active_reg_group[which(active_reg_group$Total > 1 & 
                                                          active_reg_group$Category.Total > 1),]
active_reg_multi_table = ddply(active_reg_multi,.(Category,Grouping,Category.Total),summarise,
      Exclusive=sum(Count == Total),Shared=sum(Count < Total))
active_reg_multi_table$Fraction = active_reg_multi_table$Exclusive/
  (active_reg_multi_table$Exclusive + active_reg_multi_table$Shared)
