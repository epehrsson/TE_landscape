# Load shuffled subfamily overlap for each metric

# chromHMM
print("chromHMM shuffled")
subfamily_chromHMM_shuffle = lapply(list.files(path="chromHMM/shuffled_TEs/subfamily/",pattern=".txt",full.names=TRUE)[c(1,3:10,2)],
                                    function(x) read.table(x,sep='\t'))
subfamily_chromHMM_shuffle = lapply(subfamily_chromHMM_shuffle,setNames,nm=c("subfamily","State","Sample","Length_ijk"))

subfamily_chromHMM_sample_expand = expand.grid(subfamily = rmsk_TE_subfamily$subfamily,Sample = metadata$Sample,State = chromHMM_states)
for (i in 1:10){
  subfamily_chromHMM_shuffle[[i]] = merge(subfamily_chromHMM_shuffle[[i]],subfamily_chromHMM_sample_expand,by=c("subfamily","Sample","State"),all=TRUE)
  subfamily_chromHMM_shuffle[[i]][is.na(subfamily_chromHMM_shuffle[[i]])] = 0
}
rm(subfamily_chromHMM_sample_expand)

names(subfamily_chromHMM_shuffle) = seq(1,10,1)
subfamily_chromHMM_shuffle = ldply(subfamily_chromHMM_shuffle)
colnames(subfamily_chromHMM_shuffle)[1] = "Iteration"

subfamily_chromHMM_shuffle_potential = ddply(subfamily_chromHMM_shuffle,.(Iteration,subfamily,State),summarise,Samples=sum(Length_ijk > 0))

# DNase
print("DNase shuffled")
subfamily_DNase_shuffle = lapply(list.files(path="DNase/shuffled/subfamily/",pattern="subfamily_DNase_sample_summit_",full.names=TRUE)[c(1,3:10,2)],
                                 function(x) read.table(x,sep='\t'))
subfamily_DNase_shuffle = lapply(subfamily_DNase_shuffle,setNames,nm=c("subfamily","Sample","Length_ijk"))

subfamily_DNase_expand = expand.grid(subfamily = rmsk_TE_subfamily$subfamily,Sample = as.vector(metadata[which(!is.na(metadata$DNase)),]$Sample))
for (i in 1:10){
  subfamily_DNase_shuffle[[i]] = merge(subfamily_DNase_shuffle[[i]],subfamily_DNase_expand,by=c("subfamily","Sample"),all=TRUE)
  subfamily_DNase_shuffle[[i]][is.na(subfamily_DNase_shuffle[[i]])] = 0
}
rm(subfamily_DNase_expand)

names(subfamily_DNase_shuffle) = seq(1,10,1)
subfamily_DNase_shuffle = ldply(subfamily_DNase_shuffle)
colnames(subfamily_DNase_shuffle)[1] = "Iteration"

subfamily_DNase_shuffle_potential = ddply(subfamily_DNase_shuffle,.(Iteration,subfamily),summarise,Samples=sum(Length_ijk > 0))
subfamily_DNase_shuffle_potential$State = rep("DNase",dim(subfamily_DNase_shuffle_potential)[1])

# H3K27ac
print("H3K27ac shuffled")
subfamily_H3K27ac_shuffle = lapply(list.files(path="H3K27ac/shuffled/subfamily/",pattern="subfamily_H3K27ac_sample_summit_",full.names=TRUE)[c(1,3:10,2)],
                                   function(x) read.table(x,sep='\t'))
subfamily_H3K27ac_shuffle = lapply(subfamily_H3K27ac_shuffle,setNames,nm=c("subfamily","Sample","Length_ijk"))

subfamily_H3K27ac_expand = expand.grid(subfamily = rmsk_TE_subfamily$subfamily,Sample = as.vector(metadata[which(!is.na(metadata$H3K27ac)),]$Sample))
for (i in 1:10){
  subfamily_H3K27ac_shuffle[[i]] = merge(subfamily_H3K27ac_shuffle[[i]],subfamily_H3K27ac_expand,by=c("subfamily","Sample"),all=TRUE)
  subfamily_H3K27ac_shuffle[[i]][is.na(subfamily_H3K27ac_shuffle[[i]])] = 0
}
rm(subfamily_H3K27ac_expand)

names(subfamily_H3K27ac_shuffle) = seq(1,10,1)
subfamily_H3K27ac_shuffle = ldply(subfamily_H3K27ac_shuffle)
colnames(subfamily_H3K27ac_shuffle)[1] = "Iteration"

subfamily_H3K27ac_shuffle_potential = ddply(subfamily_H3K27ac_shuffle,.(Iteration,subfamily),summarise,Samples=sum(Length_ijk > 0))
subfamily_H3K27ac_shuffle_potential$State = rep("H3K27ac",dim(subfamily_H3K27ac_shuffle_potential)[1])

# WGBS (CpGs)
print("WGBS shuffled")
subfamily_CpG_shuffle = lapply(list.files(path="WGBS/shuffled/subfamily/",pattern="subfamily_CpG_Meth_states",full.names=TRUE)[c(1,3:10,2)],
                               function(x) read.table(x,sep='\t'))
subfamily_CpG_shuffle = lapply(subfamily_CpG_shuffle,setNames,nm=c("Sample",meth_states,"subfamily"))
names(subfamily_CpG_shuffle) = seq(1,10,1)
subfamily_CpG_shuffle = ldply(subfamily_CpG_shuffle)
colnames(subfamily_CpG_shuffle)[1] = "Iteration"
subfamily_CpG_shuffle$Sample = mapvalues(subfamily_CpG_shuffle$Sample,seq(8,44,1),as.vector(metadata[which(!is.na(metadata$WGBS)),]$Sample))
subfamily_CpG_shuffle[is.na(subfamily_CpG_shuffle)] = 0
subfamily_CpG_shuffle[,meth_states] = subfamily_CpG_shuffle[,meth_states]/2
subfamily_CpG_shuffle = melt(subfamily_CpG_shuffle,id.vars=c("Iteration","Sample","subfamily"))
colnames(subfamily_CpG_shuffle)[4:5] = c("State","Length_ijk")

subfamily_CpG_shuffle_potential = ddply(subfamily_CpG_shuffle,.(Iteration,subfamily,State),summarise,Samples=sum(Length_ijk > 0))

# Combine all
print("Combine shuffled")
subfamily_state_potential_shuffle = rbind(subfamily_chromHMM_shuffle_potential,subfamily_CpG_shuffle_potential,
                                          subfamily_DNase_shuffle_potential,subfamily_H3K27ac_shuffle_potential)
subfamily_state_potential_shuffle$State = factor(subfamily_state_potential_shuffle$State,levels=states[1:21])
subfamily_state_potential_shuffle$Metric = factor(ifelse(subfamily_state_potential_shuffle$State %in% chromHMM_states,"chromHMM",
                                                         ifelse(subfamily_state_potential_shuffle$State %in% meth_states,"WGBS",
                                                                as.character(subfamily_state_potential_shuffle$State))),levels=c("chromHMM","WGBS","DNase","H3K27ac"))
subfamily_state_potential_shuffle$Sample.Proportion = as.numeric(subfamily_state_potential_shuffle$Samples/sample_counts["All",subfamily_state_potential_shuffle$Metric])
subfamily_state_potential_shuffle$Iteration = factor(subfamily_state_potential_shuffle$Iteration,levels=seq(1,10,1))

rm(list=c("subfamily_chromHMM_shuffle_potential","subfamily_CpG_shuffle_potential",
          "subfamily_DNase_shuffle_potential","subfamily_H3K27ac_shuffle_potential",
          "subfamily_chromHMM_shuffle","subfamily_CpG_shuffle","subfamily_DNase_shuffle","subfamily_H3K27ac_shuffle"))