# Proportion of chromHMM state in TEs by sample grouping
# See 5/25/2016, 9/27/2016, 2/3/2017, 3/8/2017

# Proportion of each state in all TEs in each sample
mnemonics_states_TEother_proportion = mnemonics_states_TEother/mnemonics_states

# Wilcox test for every state x Group
apply(mnemonics_states_TEother_proportion,2,function(x) wilcox_to_all(x,EID_metadata$Group))

# Mean proportion by Group for each state
apply(mnemonics_states_TEother_proportion,2,function(x) tapply(x,EID_metadata$Group,mean))

mnemonics_states_TEother_proportion_long = melt(as.matrix(mnemonics_states_TEother_proportion))
colnames(mnemonics_states_TEother_proportion_long) = c("Sample","State","Proportion")
mnemonics_states_TEother_proportion_long$Type = rep(EID_metadata$Type,16)
mnemonics_states_TEother_proportion_long$Anatomy = rep(EID_metadata$Anatomy,16)
mnemonics_states_TEother_proportion_long$Group = rep(EID_metadata$Group,16)
mnemonics_states_TEother_proportion_long$Age = rep(EID_metadata$Age,16)
mnemonics_states_TEother_proportion_long$Cancer = rep(EID_metadata$Cancer,16)
mnemonics_states_TEother_proportion_long$Germline = rep(EID_metadata$Germline,16)

# Min, max, median of proportion of each state in all TEs in each sample
mnemonics_states_TEother_proportion_range = rbind(apply(mnemonics_states_TEother_proportion,2,min),apply(mnemonics_states_TEother_proportion,2,max),apply(mnemonics_states_TEother_proportion,2,median))
rownames(mnemonics_states_TEother_proportion_range) = c("Min","Max","Median")

# Kruskal-Wallis with Bonferroni for signficantly different state proportion in all TEs between sample groups
mnemonics_states_TEother_proportion_kruskal = rbind(p.adjust(apply(mnemonics_states_TEother_proportion,2,function(x) unlist(kruskal.test(x ~ EID_metadata$Group))["p.value"]),method="bonf"),p.adjust(apply(mnemonics_states_TEother_proportion,2,function(x) unlist(kruskal.test(x ~ EID_metadata$Anatomy))["p.value"]),method="bonf"),p.adjust(apply(mnemonics_states_TEother_proportion,2,function(x) unlist(kruskal.test(x ~ EID_metadata$Type))["p.value"]),method="bonf"),p.adjust(apply(mnemonics_states_TEother_proportion,2,function(x) unlist(kruskal.test(x ~ EID_metadata$Cancer))["p.value"]),method="bonf"),p.adjust(apply(mnemonics_states_TEother_proportion,2,function(x) unlist(kruskal.test(x ~ EID_metadata$Age))["p.value"]),method="bonf"),p.adjust(apply(mnemonics_states_TEother_proportion,2,function(x) unlist(kruskal.test(x ~ EID_metadata$Germline))["p.value"]),method="bonf"))
rownames(mnemonics_states_TEother_proportion_kruskal) = c("Group","Anatomy","Type","Cancer","Age","Germline")
