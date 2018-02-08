test = aggregate(data=subfamily_state_sample_combined,Members~subfamily+State,function(x) sum(x > 1))
samples = lapply(c(chromHMM_states,meth_states,"DNase","H3K27ac"),
                 function(x) as.vector(test[which(test$State == x & test$Members > 1),]$subfamily))
names(samples) = c(chromHMM_states,meth_states,"DNase","H3K27ac")
lapply(samples$'5_TxWk'[c(223:299,308:762,767:961)],function(x) investigate_candidate_indv(subfamily=x,state="5_TxWk",metric="chromHMM",print_fig=TRUE))
lapply(chromHMM_states[1:7],function(y) lapply(samples[[y]],function(x) investigate_candidate_indv(subfamily=x,state=y,metric="chromHMM",print_fig=TRUE)))
