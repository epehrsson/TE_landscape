# Total bases/CpGs in state by sample, genome and TEs
# Requires proportion.R

# chromHMM
chromHMM_state_proportion = rbind(melt(as.matrix(mnemonics_states_genome[2:16,1:127])),melt(as.matrix(mnemonics_states_TE[2:16,1:127])))
colnames(chromHMM_state_proportion) = c("State","Sample","Bases")
chromHMM_state_proportion$Cohort = rep(c("Genome","TE"),each=1905)
  
# WGBS
WGBS_state_proportion = rbind(melt(as.matrix(all_CpG_meth)),melt(as.matrix(TE_CpG_meth)))
colnames(WGBS_state_proportion) = c("Sample","State","CpGs")
WGBS_state_proportion$Cohort = rep(c("Genome","TE"),each=148)

# DNase
DNase_stats_long = melt(DNase_stats[,c(1,4:5)],id.var=("Sample"))
colnames(DNase_stats_long)[2:3] = c("Cohort","Bases")
DNase_stats_long$State = rep("DNase",dim(DNase_stats_long)[1])

# H3K27ac
H3K27ac_stats_long = melt(H3K27ac_stats[,c(1,4:5)],id.var=("Sample"))
colnames(H3K27ac_stats_long)[2:3] = c("Cohort","Bases")
H3K27ac_stats_long$State = rep("H3K27ac",dim(H3K27ac_stats_long)[1])

# Combined
test = WGBS_state_proportion
colnames(test)[3] = "Bases"
all_state_proportion = rbind(chromHMM_state_proportion,test,DNase_stats_long,H3K27ac_stats_long)
all_state_proportion$Cohort = gsub("Total_width_in_TE","TE",all_state_proportion$Cohort)
all_state_proportion$Cohort = gsub("Total_width","Genome",all_state_proportion$Cohort)
all_state_proportion$Mark = factor(c(rep("chromHMM",3810),rep("WGBS",296),rep("DNase",106),rep("H3K27ac",196)),levels=c("chromHMM","WGBS","DNase","H3K27ac"))
rm(test)

# Add metadata
all_state_proportion = merge(all_state_proportion,metadata[,c(1,4:9)],by="Sample")
