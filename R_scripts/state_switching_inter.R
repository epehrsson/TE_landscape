# Probability that a TE ever in a particular state will be in another state in a different sample
# See 7/7/2016, 2/6/2017, 2/9/2017

#Format TEs
ss_inter_TE = read.table("all_chromHMM_TE_state_inter.txt",sep='\t',header=TRUE,row.names=1)
ss_inter_TE[lower.tri(ss_inter_TE)] = t(ss_inter_TE)[lower.tri(ss_inter_TE)]

#Format other TEs
ss_inter_other = read.table("state_switching/all_chromHMM_other_state_inter.txt",sep='\t',header=TRUE,row.names=1)
ss_inter_other[lower.tri(ss_inter_other)] = t(ss_inter_other)[lower.tri(ss_inter_other)]

#Combine
state_switching_inter = ss_inter_TE + ss_inter_other

#Normalize
for (i in 1:15){
  state_switching_inter[i,] = state_switching_inter[i,]/apply(potential_TEother_state_dist[,2:16],2,function(x) sum(x[2:128]))[i]
  # Was state_switching_inter[i,] = state_switching_inter[i,]/apply(potential_TE_state_ratchet[,2:16],2,function(x) sum(x[2:128]))[i]
}

rm(ss_inter_TE)
rm(ss_inter_other)