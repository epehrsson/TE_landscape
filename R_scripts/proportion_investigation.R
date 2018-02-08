# chromHMM
# Combine genome/TE, changing colnames, adding column of cohort
chromHMM_state_proportion = rbind(melt(as.matrix(mnemonics_states_genome[2:16,1:127])),melt(as.matrix(mnemonics_states_TE[2:16,1:127])))
colnames(chromHMM_state_proportion) = c("State","Sample","Bases")
chromHMM_state_proportion$Cohort = rep(c("Genome","TE"),each=1905)

# Add class
test = class_chromHMM[,1:4]
colnames(test)[c(2,4)] = c("Cohort","Bases")
chromHMM_state_proportion = rbind(chromHMM_state_proportion,test)
chromHMM_state_proportion$Cohort = factor(chromHMM_state_proportion$Cohort,levels=c("Genome","TE","DNA","LINE","LTR","SINE","SVA","Other"))

# Merge with metadata
chromHMM_state_proportion = merge(chromHMM_state_proportion,metadata[,c(1,4)],by="Sample")
  
# Plot genome and TE bases in state per sample
lapply(chromHMM_states,function(x) ggplot(chromHMM_state_proportion[which(chromHMM_state_proportion$Cohort %in% c("Genome","TE") & chromHMM_state_proportion$State == x),],aes(x=reorder(Sample,Bases,max),y=log10(Bases),color=Group,group=Cohort,shape=Cohort)) + 
         geom_point(size=2) + geom_line() + theme(axis.text.x = element_text(angle=90),axis.title.x = element_blank()) + scale_color_manual(values=group_colors))

# Compare TE proportion to genome, non-TE, and TE bases
chromHMM_state_proportion_matrix = dcast(chromHMM_state_proportion,State+Sample+Group~Cohort,value.var="Bases")
ddply(chromHMM_state_proportion_matrix,.(State),summarise,Genome_corr = unlist(cor.test((TE/Genome),Genome))["estimate.cor"],Genome_pvalue = unlist(cor.test((TE/Genome),Genome))["p.value"],
      TE_corr = unlist(cor.test((TE/Genome),TE))["estimate.cor"],TE_pvalue = unlist(cor.test((TE/Genome),TE))["p.value"],
      NonTE_corr = unlist(cor.test((TE/Genome),(Genome-TE)))["estimate.cor"],NonTE_pvalue = unlist(cor.test((TE/Genome),(Genome-TE)))["p.value"])

lapply(chromHMM_states,function(x) ggplot(chromHMM_state_proportion_matrix[which(chromHMM_state_proportion_matrix$State == x),],aes(x=reorder(Sample,TE/Genome,max),y=(TE/(Genome-TE)),color=Group)) + geom_point(shape=16) + scale_color_manual(values=group_colors) + theme(axis.text.x = element_text(angle=90),axis.title.x = element_blank()))

# WGBS
WGBS_state_proportion = rbind(melt(as.matrix(all_CpG_meth)),melt(as.matrix(TE_CpG_meth)))
colnames(WGBS_state_proportion) = c("Sample","State","CpGs")
WGBS_state_proportion$Cohort = rep(c("Genome","TE"),each=148)

# Merge with metadata
WGBS_state_proportion = merge(WGBS_state_proportion,metadata[,c(1,4)],by="Sample")

# Plot genome and TE CpGs in state per sample
lapply(meth_states,function(x) ggplot(WGBS_state_proportion[which(WGBS_state_proportion$State == x),],aes(x=reorder(Sample,CpGs,max),y=CpGs,color=Group,group=Cohort,shape=Cohort)) + 
         geom_point(size=2) + geom_line() + theme(axis.text.x = element_text(angle=90),axis.title.x = element_blank()) + scale_color_manual(values=group_colors))

WGBS_state_proportion_matrix = dcast(WGBS_state_proportion,State+Sample+Group~Cohort,value.var="CpGs")
ddply(WGBS_state_proportion_matrix,.(State),summarise,Genome_corr = unlist(cor.test((TE/Genome),Genome))["estimate.cor"],Genome_pvalue = unlist(cor.test((TE/Genome),Genome))["p.value"],
      TE_corr = unlist(cor.test((TE/Genome),TE))["estimate.cor"],TE_pvalue = unlist(cor.test((TE/Genome),TE))["p.value"],
      NonTE_corr = unlist(cor.test((TE/Genome),(Genome-TE)))["estimate.cor"],NonTE_pvalue = unlist(cor.test((TE/Genome),(Genome-TE)))["p.value"])

# DNase
DNase_stats_long = melt(DNase_stats[,c(1,4:5)],id.var=("Sample"))
colnames(DNase_stats_long)[2:3] = c("Cohort","Bases")
DNase_stats_long = merge(DNase_stats_long,metadata[,c(1,4:5)],by="Sample")

ggplot(DNase_stats_long,aes(x=reorder(Sample,Bases,max),y=Bases,color=Group,group=Cohort,shape=Cohort)) + geom_point(size=2) + geom_line() + theme(axis.text.x = element_text(angle=90),axis.title.x = element_blank()) + scale_color_manual(values=group_colors)

cor.test(DNase_stats$Total_width,(DNase_stats$Total_width_in_TE/DNase_stats$Total_width))
cor.test(DNase_stats$Total_width-DNase_stats$Total_width_in_TE,(DNase_stats$Total_width_in_TE/DNase_stats$Total_width))
cor.test(DNase_stats$Total_width_in_TE,(DNase_stats$Total_width_in_TE/DNase_stats$Total_width))

# H3K27ac
H3K27ac_stats_long = melt(H3K27ac_stats[,c(1,4:5)],id.var=("Sample"))
colnames(H3K27ac_stats_long)[2:3] = c("Cohort","Bases")
H3K27ac_stats_long = merge(H3K27ac_stats_long,metadata[,c(1,4:5)],by="Sample")

ggplot(H3K27ac_stats_long,aes(x=reorder(Sample,Bases,max),y=Bases,color=Group,group=Cohort,shape=Cohort)) + geom_point(size=2) + geom_line() + theme(axis.text.x = element_text(angle=90),axis.title.x = element_blank()) + scale_color_manual(values=group_colors)

cor.test(H3K27ac_stats$Total_width,(H3K27ac_stats$Total_width_in_TE/H3K27ac_stats$Total_width))
cor.test(H3K27ac_stats$Total_width-H3K27ac_stats$Total_width_in_TE,(H3K27ac_stats$Total_width_in_TE/H3K27ac_stats$Total_width))
cor.test(H3K27ac_stats$Total_width_in_TE,(H3K27ac_stats$Total_width_in_TE/H3K27ac_stats$Total_width))
