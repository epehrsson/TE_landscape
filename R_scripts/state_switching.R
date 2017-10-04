# State switching
# See 7/5/2016, 7/7/2016, 2/6/2017, 2/9/2017

# Matrix of states annotated on the same TE in the same sample (intra)
# Normalized by number of instances (TE x sample) each state is found, alone or with another state

source("R_scripts/potential_class.R")

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
rm(ss_intra_TE)
rm(ss_intra_other)

# For Chi-sq tests
state_switching_intra_raw = state_switching_intra

#Normalize matrix
for (i in 1:15){
  state_switching_intra[i,] = state_switching_intra[i,]/state_switching_intra[i,i]
}

# chromHMM inter state switching
# Format TEs
ss_inter_TE = read.table("chromHMM/state_switching/all_chromHMM_TE_state_inter.txt",sep='\t',header=TRUE,row.names=1)
ss_inter_TE[lower.tri(ss_inter_TE)] = t(ss_inter_TE)[lower.tri(ss_inter_TE)]

# Format other TEs
ss_inter_other = read.table("chromHMM/state_switching/all_chromHMM_other_state_inter.txt",sep='\t',header=TRUE,row.names=1)
ss_inter_other[lower.tri(ss_inter_other)] = t(ss_inter_other)[lower.tri(ss_inter_other)]

# Combine
state_switching_inter = ss_inter_TE + ss_inter_other
rm(ss_inter_TE)
rm(ss_inter_other)

# For Chi-sq tests
state_switching_inter_raw = state_switching_inter

# Normalize
for (i in 1:15){
  state_switching_inter[i,] = state_switching_inter[i,]/apply(chromHMM_TE_state_dist[,2:16],2,function(x) sum(x[2:128]))[i]
}

# chromHMM, by class
ss_inter_class = lapply(list.files(path="chromHMM/state_switching/class",full.names=TRUE),function(x) read.table(x,sep='\t',header=TRUE,row.names=1))
names(ss_inter_class) = c("DNA","LINE","LTR","Other","SINE","SVA")

for (j in 1:6){
 ss_inter_class[[j]][lower.tri(ss_inter_class[[j]])] = t(ss_inter_class[[j]])[lower.tri(ss_inter_class[[j]])]
 for (i in 1:15){
   ss_inter_class[[j]][i,] = ss_inter_class[[j]][i,]/apply(chromHMM_TE_state_class[which(chromHMM_TE_state_class$class_update == names(ss_inter_class)[j]),3:17],2,function(x) sum(x[2:128]))[i]
 }
}

ss_inter_class = ldply(ss_inter_class)
colnames(ss_inter_class)[1] = c("Class")
ss_inter_class$State = rep(chromHMM_states,6)
ss_inter_class$Class = factor(ss_inter_class$Class,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))

# WGBS inter state switching
# Format TEs
ss_inter_meth = read.table("WGBS/TE_WGBS_state_inter.txt",sep='\t',header=TRUE,row.names=1)
ss_inter_meth[lower.tri(ss_inter_meth)] = t(ss_inter_meth)[lower.tri(ss_inter_meth)]

# Reorder columns and rows
ss_inter_meth = ss_inter_meth[c(2,3,1,4),c(2,3,1,4)]

# For Chi-sq tests
ss_inter_meth_raw = ss_inter_meth

# Normalize
for (i in 1:4){
  ss_inter_meth[i,] = ss_inter_meth[i,]/apply(TE_meth_average_category[,c(2,4,3,5)],2,function(x) sum(x[2:38]))[i]
}

# WGBS, by class
ss_inter_meth_class = lapply(list.files(path="WGBS/class",pattern="WGBS_state_sorted.txt_inter.txt",full.names=TRUE),function(x) read.table(x,sep='\t',header=TRUE,row.names=1))
names(ss_inter_meth_class) = c("DNA","LINE","LTR","Other","SINE","SVA")

for (j in 1:6){
  ss_inter_meth_class[[j]][lower.tri(ss_inter_meth_class[[j]])] = t(ss_inter_meth_class[[j]])[lower.tri(ss_inter_meth_class[[j]])]
  ss_inter_meth_class[[j]] = ss_inter_meth_class[[j]][c(2,3,1,4),c(2,3,1,4)]
  for (i in 1:4){
    ss_inter_meth_class[[j]][i,] = ss_inter_meth_class[[j]][i,]/apply(TE_meth_average_class[which(TE_meth_average_class$class_update == names(ss_inter_meth_class)[j]),c(3,5,4,6)],2,function(x) sum(x[2:38]))[i]
  }
}
ss_inter_meth_class = ldply(ss_inter_meth_class)
colnames(ss_inter_meth_class)[1] = c("Class")
ss_inter_meth_class$State = rep(meth_states,6)
ss_inter_meth_class$Class = factor(ss_inter_meth_class$Class,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
