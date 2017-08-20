# Probability that a TE ever in a particular state will be in another state in a different sample
# See 7/5/2016, 7/7/2016, 2/6/2017, 2/9/2017

#Format TEs
ss_inter_TE = read.table("chromHMM/state_switching/all_chromHMM_TE_state_inter.txt",sep='\t',header=TRUE,row.names=1)
ss_inter_TE[lower.tri(ss_inter_TE)] = t(ss_inter_TE)[lower.tri(ss_inter_TE)]

#Format other TEs
ss_inter_other = read.table("chromHMM/state_switching/all_chromHMM_other_state_inter.txt",sep='\t',header=TRUE,row.names=1)
ss_inter_other[lower.tri(ss_inter_other)] = t(ss_inter_other)[lower.tri(ss_inter_other)]

#Combine
state_switching_inter = ss_inter_TE + ss_inter_other

#Normalize
for (i in 1:15){
  state_switching_inter[i,] = state_switching_inter[i,]/apply(chromHMM_TE_state_dist[,2:16],2,function(x) sum(x[2:128]))[i]
}

rm(ss_inter_TE)
rm(ss_inter_other)

# Matrix of states annotated on the same TE in the same sample
# Normalized by number of instances (TE x sample) each state is found, alone or with another state

# Format TEs 
ss_intra_TE = read.table("chromHMM/state_switching/all_chromHMM_TE_state_intra.txt",sep='\t',header=TRUE,row.names=1)
ss_intra_TE[1:9,10:15] = t(ss_intra_TE[10:15,1:9])
ss_intra_TE[lower.tri(ss_intra_TE)] = t(ss_intra_TE)[lower.tri(ss_intra_TE)]

# Format other TEs
ss_intra_other = read.table("chromHMM/state_switching/all_chromHMM_other_state_intra.txt",sep='\t',header=TRUE,row.names=1)
ss_intra_other[1:9,10:15] = t(ss_intra_other[10:15,1:9])
ss_intra_other[lower.tri(ss_intra_other)] = t(ss_intra_other)[lower.tri(ss_intra_other)]

# Combine
state_switching_intra = ss_intra_TE + ss_intra_other

#Normalize matrix
for (i in 1:15){
  state_switching_intra[i,] = state_switching_intra[i,]/state_switching_intra[i,i]
}

rm(ss_intra_TE)
rm(ss_intra_other)