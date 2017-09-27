# Proportion of each state in TEs by sample grouping

source("R_scripts/proportion.R")
source("R_scripts/proportion_class.R")

# chromHMM
# Proportion of each chromHMM state in TEs by sample
mnemonics_states_TE_proportion = melt(as.matrix(mnemonics_states_TE[2:16,1:127]/mnemonics_states_genome[2:16,1:127]))
colnames(mnemonics_states_TE_proportion) = c("State","Sample","Proportion")

# WGBS
# Proportion of hypomethylated CpGs in TEs by sample
WGBS_proportion = melt(as.matrix(TE_CpG_meth/all_CpG_meth))
colnames(WGBS_proportion) = c("Sample","State","Proportion")

# DNase
# Proportion of Dnase peaks overlapping TEs by sample classification
DNase_proportion = as.data.frame(cbind(as.vector(DNase_stats$Sample),DNase_stats$Total_width_in_TE/DNase_stats$Total_width))
colnames(DNase_proportion) = c("Sample","Proportion")
DNase_proportion$State = rep("DNase",53)

# H3K27ac
# Proportion of H3K27ac peaks overlapping TEs by sample classification
H3K27ac_proportion = as.data.frame(cbind(as.vector(H3K27ac_stats$Sample),H3K27ac_stats$Total_width_in_TE/H3K27ac_stats$Total_width))
colnames(H3K27ac_proportion) = c("Sample","Proportion")
H3K27ac_proportion$State = rep("H3K27ac",98)

# By class
# chromHMM
class_sample_chromHMM = class_chromHMM[,1:4]
class_sample_chromHMM$Proportion = apply(class_sample_chromHMM,1,function(x) as.numeric(x[4])/mnemonics_states_genome[x[3],x[1]])

# WGBS
class_sample_WGBS = class_CpG_meth[,1:4]
class_sample_WGBS$Proportion = apply(class_sample_WGBS,1,function(x) as.numeric(x[4])/all_CpG_meth[x[2],x[3]])

# DNase
class_sample_DNase = TE_DNase_class[,c(1:3,5)]
class_sample_DNase$Proportion = apply(class_sample_DNase,1,function(x) as.numeric(x[3])/DNase_stats[match(x[2],DNase_stats$Sample),]$Total_width)

# H3K27ac
class_sample_H3K27ac = TE_H3K27ac_class[,c(1:3,5)]
class_sample_H3K27ac$Proportion = apply(class_sample_H3K27ac,1,function(x) as.numeric(x[3])/H3K27ac_stats[match(x[2],H3K27ac_stats$Sample),]$Total_width)
