# For 10 iterations of shuffled TEs, chromHMM and methylation states, creates matrices with:
# 1) the number of unique states each TE is annotated with across samples
# 2) the average number of samples a TE is annotated with each state, for those ever annotated with each state
# 3) the median number of states each TE is annotated with across samples, for TEs ever annotated with each state

## shuffled_chromHMM_states - number of unique states each TE is annotated with across samples, chromHMM
## shuffled_WGBS_states - number of unique states each TE is annotated with across samples, methylation
## shuffled_chromHMM_inter - average number of samples a TE is annotated with each state, for those ever annotated with each state, chromHMM
## shuffled_meth_inter - average number of samples a TE is annotated with each state, for those ever annotated with each state, methylation
## shuffled_chromHMM_dynamics - median number of states each TE is annotated with across all samples, for TEs ever in each state, chromHMM
## shuffled_meth_dynamics - median number of states each TE is annotated with across all samples, for TEs ever in each state, methylation

# Load dataframes with:
# chromHMM: the number of samples in which each shuffled TE is annotated with each state
# and the number of unique states the TE is annotated with across all samples
# WGBS: the average methylation per shuffled TE for each sample, for only TEs that overlap at least one CpG,
# the number of CpGs per TE, the number of samples in which the TE is annotated with each methylation state, 
# and the total number of methylation states with which the TE is annotated across samples

load("R_datasets/shuffled_chromHMM.RData")
load("R_datasets/shuffled_WGBS.RData")


# Dataframes with the number of unique states each TE is annotated with across samples, plus class annotation

# chromHMM states
shuffled_chromHMM_states = ldply(shuffled_chromHMM_potential,function(x) x[,c("class","States")])
## Add iteration
shuffled_chromHMM_states$Iteration = factor(rep(seq(1,10,1),each=NUM_TE))

# Methylation states
## Add iteration
shuffled_WGBS_average = lapply(seq(1,10,1),function(y) transform(shuffled_WGBS_average[[y]],Iteration = rep(y,unlist(dim(shuffled_WGBS_average[[y]])[1]))))
shuffled_WGBS_states = ldply(shuffled_WGBS_average,function(x) x[,c("class","States","Iteration")])
shuffled_WGBS_states$Iteration = as.factor(shuffled_WGBS_states$Iteration)


# Dataframes with the average number of samples a TE is annotated with each state, for those ever annotated with each state

# Number of TEs annotated with each chromHMM state in at least one sample
shuffled_chromHMM_ever = ldply(shuffled_chromHMM_potential,function(x) apply(x[,8:22],2,function(y) sum(y > 0)))
colnames(shuffled_chromHMM_ever) = chromHMM_states
shuffled_chromHMM_ever$Iteration = seq(1,10,1)
shuffled_chromHMM_ever = melt(shuffled_chromHMM_ever,id.vars="Iteration")
colnames(shuffled_chromHMM_ever)[2:3] = c("State","Count")

# Load matrices of the total number of samples in which TEs are annotated with each chromHMM state (column)
# across all TEs annotated with each chromHMM state in any sample (row)
shuffled_chromHMM_inter = lapply(seq(1,10,1),function(x) read.table(paste("chromHMM/shuffled_TEs/rmsk_TE_shuffle_",x,"_chromHMM_inter.txt",sep=""),sep='\t',header=TRUE,row.names=1))

## Reformat
shuffled_chromHMM_inter = ldply(shuffled_chromHMM_inter,function(x) melt(as.matrix(x)))
colnames(shuffled_chromHMM_inter) = c("State1","State2","Count")
shuffled_chromHMM_inter$Iteration = rep(seq(1,10,1),each=225)

## Normalize each row by the number of TEs ever in the state
shuffled_chromHMM_inter = merge(shuffled_chromHMM_inter,shuffled_chromHMM_ever,by.x=c("Iteration","State1"),by.y=c("Iteration","State"))
shuffled_chromHMM_inter$Mean_samples = shuffled_chromHMM_inter$Count.x/shuffled_chromHMM_inter$Count.y

# Load matrices of the total number of samples in which TEs are annotated with each methylation state (column)
# across all TEs annotated with each methylation state in any sample (row)
# Plus the number of TEs ever in each state (Total)
shuffled_meth_inter = lapply(seq(1,10,1),function(x) read.table(paste("WGBS/shuffled/TE_shuffle_",x,"_inter.txt",sep=""),sep='\t',header=TRUE,row.names=1))

## Reformat
shuffled_meth_inter = lapply(shuffled_meth_inter,function(x) transform(x,State1=rownames(x)))
shuffled_meth_inter = ldply(shuffled_meth_inter,function(x) melt(x,id.vars=c("State1","Total")))
colnames(shuffled_meth_inter)[3:4] = c("State2","Count")
shuffled_meth_inter$Iteration = rep(seq(1,10,1),each=16)

## Normalize each row by the number of TEs ever in the state
shuffled_meth_inter$Mean_samples = shuffled_meth_inter$Count/shuffled_meth_inter$Total


# Dataframes with the median number of states each TE is annotated with across all samples, for TEs ever in each state

# chromHMM states
shuffled_chromHMM_dynamics = lapply(shuffled_chromHMM_potential,function(y) apply(y[,8:22],2,function(x) median(y[which(x > 0),]$States)))
shuffled_chromHMM_dynamics = ldply(shuffled_chromHMM_dynamics)

# Methylation states
shuffled_meth_dynamics = lapply(shuffled_WGBS_average,function(y) apply(y[,meth_states],2,function(x) median(y[which(x > 0),]$States)))
shuffled_meth_dynamics = ldply(shuffled_meth_dynamics)

rm(list=c("shuffled_chromHMM_ever","shuffled_WGBS_average"))