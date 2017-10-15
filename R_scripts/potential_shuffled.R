# Potential for shuffled matrices
 
library(reshape2)

load("R_datasets/shuffled_chromHMM.RData")
load("R_datasets/shuffled_WGBS.RData")
load("R_datasets/shuffled_DNase.RData")
load("R_datasets/shuffled_H3K27ac.RData")

# chromHMM
shuffled_chromHMM_potential_stats = lapply(shuffled_chromHMM_potential,function(x) potential_stats(sample_distribution(x,c(8:22),127),15,127))
shuffled_chromHMM_potential_stats = lapply(shuffled_chromHMM_potential_stats,function(x) transform(x,State = factor(chromHMM_states,levels=chromHMM_states)))
shuffled_chromHMM_potential_stats = ldply(shuffled_chromHMM_potential_stats)
shuffled_chromHMM_potential_stats$Iteration = rep(seq(1,10,1),each=15)

# WGBS
shuffled_WGBS_average_stats = lapply(shuffled_WGBS_average,function(x) potential_stats(sample_distribution(x,c(46:49),37),4,37))
shuffled_WGBS_average_stats = lapply(shuffled_WGBS_average_stats,function(x) transform(x,State = factor(rownames(x),levels=meth_states)))
shuffled_WGBS_average_stats = ldply(shuffled_WGBS_average_stats)
shuffled_WGBS_average_stats$Iteration = rep(seq(1,10,1),each=4)

# WGBS, no IMR90
shuffled_WGBS_average_noIMR90_stats = lapply(shuffled_WGBS_average,function(x) potential_stats(sample_distribution(x,c(50:53),36),4,36))
shuffled_WGBS_average_noIMR90_stats = lapply(shuffled_WGBS_average_noIMR90_stats,function(x) transform(x,State = factor(c("Hypomethylated","Hypermethylated","Intermediate","Missing"),levels=meth_states)))
shuffled_WGBS_average_noIMR90_stats = ldply(shuffled_WGBS_average_noIMR90_stats)
shuffled_WGBS_average_noIMR90_stats$Iteration = rep(seq(1,10,1),each=4)

# DNase
shuffled_DNase_potential_dist = lapply(shuffled_DNase_potential,function(x) sample_distribution(x,61,53))
shuffled_DNase_potential_dist = lapply(shuffled_DNase_potential_dist,setNames,nm=c("Samples","DNase"))
for (i in 1:10){
  shuffled_DNase_potential_dist[[i]][1,2] = 4430788-dim(shuffled_DNase_potential[[i]])[1]
}
shuffled_DNase_potential_stats = lapply(shuffled_DNase_potential_dist,function(x) potential_stats(x,1,53))
shuffled_DNase_potential_stats = lapply(shuffled_DNase_potential_stats,function(x) transform(x,State = "DNase"))
shuffled_DNase_potential_stats = ldply(shuffled_DNase_potential_stats)
shuffled_DNase_potential_stats$Iteration = seq(1,10,1)

# H3K27ac 
shuffled_H3K27ac_potential_dist = lapply(shuffled_H3K27ac_potential,function(x) sample_distribution(x,106,98))
shuffled_H3K27ac_potential_dist = lapply(shuffled_H3K27ac_potential_dist,setNames,nm=c("Samples","H3K27ac"))
for (i in 1:10){
  shuffled_H3K27ac_potential_dist[[i]][1,2] = 4430788-dim(shuffled_H3K27ac_potential[[i]])[1]
}
shuffled_H3K27ac_potential_stats = lapply(shuffled_H3K27ac_potential_dist,function(x) potential_stats(x,1,98))
shuffled_H3K27ac_potential_stats = lapply(shuffled_H3K27ac_potential_stats,function(x) transform(x,State = "H3K27ac"))
shuffled_H3K27ac_potential_stats = ldply(shuffled_H3K27ac_potential_stats)
shuffled_H3K27ac_potential_stats$Iteration = seq(1,10,1)

# Number of shuffled TEs in state by sample

# chromHMM
state_sample_count_shuffled = lapply(list.files(path="chromHMM/shuffled_TEs/",pattern="state_sample_count_",full.names = TRUE),function(x) read.table(x,sep='\t'))
state_sample_count_shuffled = lapply(state_sample_count_shuffled, setNames, nm =c("Sample","State","Count"))
state_sample_count_shuffled = lapply(state_sample_count_shuffled,function(y) transform(y,Proportion = apply(y,1,function(x) {ifelse(metadata[match(x[1],metadata$Sample),]$chrY == "Yes",as.numeric(x[3])/4430788,as.numeric(x[3])/(4430788-31580))})))
state_sample_count_shuffled = lapply(state_sample_count_shuffled,function(y) transform(y,State = factor(y$State,levels=chromHMM_states)))
state_sample_count_shuffled = ldply(state_sample_count_shuffled)
state_sample_count_shuffled$Iteration = rep(seq(1,10,1),each=1905)

# WGBS
shuffled_WGBS_average_count = lapply(shuffled_WGBS_average,function(y) as.data.frame(cbind(apply(y[,8:44],2,function(x) sum(na.omit(x) < 0.3)/length(x)),apply(y[,8:44],2,function(x)  sum(na.omit(x) <= 0.7 & na.omit(x) >= 0.3)/length(x)),apply(y[,8:44],2,function(x) sum(na.omit(x) > 0.7)/length(x)),apply(y[,8:44],2,function(x) sum(is.na(x))/length(x)))))
shuffled_WGBS_average_count = lapply(shuffled_WGBS_average_count, setNames, nm=c("Hypomethylated","Intermediate","Hypermethylated","Missing"))
shuffled_WGBS_average_count = lapply(shuffled_WGBS_average_count,function(y) melt(as.matrix(y)))
shuffled_WGBS_average_count = lapply(shuffled_WGBS_average_count, setNames, nm=c("Sample","State","Proportion"))
shuffled_WGBS_average_count = ldply(shuffled_WGBS_average_count)
shuffled_WGBS_average_count$Iteration = rep(seq(1,10,1),each=148)

# DNase
shuffled_DNase_potential_count = lapply(shuffled_DNase_potential,function(y) as.data.frame(apply(y[,8:60],2,function(x) sum(x > 0))))
shuffled_DNase_potential_count = lapply(shuffled_DNase_potential_count, setNames, nm="Count")
shuffled_DNase_potential_count = lapply(shuffled_DNase_potential_count,function(y) transform(y,Sample = rownames(y),State = rep("DNase",53)))
shuffled_DNase_potential_count = lapply(shuffled_DNase_potential_count,function(y) transform(y,Proportion = apply(y,1,function(x) {ifelse(metadata[match(x[2],metadata$Sample),]$chrY == "Yes",as.numeric(x[1])/4430788,as.numeric(x[1])/(4430788-31580))})))
shuffled_DNase_potential_count = ldply(shuffled_DNase_potential_count)
shuffled_DNase_potential_count$Iteration = rep(seq(1,10,1),each=53)

# H3K27ac
shuffled_H3K27ac_potential_count = lapply(shuffled_H3K27ac_potential,function(y) as.data.frame(apply(y[,8:105],2,function(x) sum(x > 0))))
shuffled_H3K27ac_potential_count = lapply(shuffled_H3K27ac_potential_count, setNames, nm="Count")
shuffled_H3K27ac_potential_count = lapply(shuffled_H3K27ac_potential_count,function(y) transform(y,Sample = rownames(y),State = rep("H3K27ac",98)))
shuffled_H3K27ac_potential_count = lapply(shuffled_H3K27ac_potential_count,function(y) transform(y,Proportion = apply(y,1,function(x) {ifelse(metadata[match(x[2],metadata$Sample),]$chrY == "Yes",as.numeric(x[1])/4430788,as.numeric(x[1])/(4430788-31580))})))
shuffled_H3K27ac_potential_count = ldply(shuffled_H3K27ac_potential_count)
shuffled_H3K27ac_potential_count$Iteration = rep(seq(1,10,1),each=98)
