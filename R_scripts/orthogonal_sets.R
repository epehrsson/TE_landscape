# Restricted TEs
restricted = lapply(states[c(1:3,6:7)],function(x) merge(chromHMM_TE_state[which(chromHMM_TE_state[[x]] > 0 & chromHMM_TE_state[[x]] <= 5),TE_coordinates],rmsk_TE,by=TE_coordinates))
names(restricted) = states[c(1:3,6:7)]

lapply(states[c(1:3,6:7)],function(x) write.table(restricted[[x]][,TE_coordinates[c(1:4,6,7,5)]],file=paste("orthogonal/restricted_",x,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t'))

# Positive control, shuffled TEs
load("R_datasets/shuffled_TEs.RData")
load("R_datasets/shuffled_chromHMM.RData")
shuffled_pos = lapply(seq(1,10,1),
                      function(y) {a = lapply(states[c(1:3,6:7)],function(x) merge(shuffled_chromHMM_potential[[y]][which(shuffled_chromHMM_potential[[y]][[paste("X",x,sep="")]] > 50),TE_coordinates],
                                                                                   rmsk_TE_shuffled[[y]],by=TE_coordinates));
                        names(a) = states[c(1:3,6:7)]; return(a)})
rm(list=c("shuffled_chromHMM_potential","rmsk_TE_shuffled"))

lapply(seq(1,10,1),function(y) lapply(states[c(1:3,6:7)],function(x) write.table(shuffled_pos[[y]][[x]][,TE_coordinates[c(1:4,6,7,5)]],
                                                                                 file=paste("orthogonal/shuffle_",y,"_constitutive_",x,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')))

# Positive control, TEs
TE_pos = lapply(states[c(1:3,6:7)],function(x) merge(chromHMM_TE_state[which(chromHMM_TE_state[[x]] > 50),c(TE_coordinates)],rmsk_TE,by=TE_coordinates))
names(TE_pos) = states[c(1:3,6:7)]

lapply(states[c(1:3,6:7)],function(x) write.table(TE_pos[[x]][,TE_coordinates[c(1:4,6,7,5)]],file=paste("orthogonal/constitutive_",x,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t'))

# Negative control
TE_neg = merge(chromHMM_TE_state[which(apply(chromHMM_TE_state[,states[c(1:3,6:7)]],1,function(x) sum(x) == 0)),c(TE_coordinates)],rmsk_TE,by=TE_coordinates)

write.table(TE_neg[,TE_coordinates[c(1:4,6,7,5)]],file="orthogonal/not_active.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')

# Combine profiles
profile_set = function(matrix_list){
  profiled = ldply(matrix_list,function(x) c(dim(x)[1],
                                             median(x$stop-x$start),
                                             table(x$class)/dim(x)[1],
                                             apply(x[,cohorts[c(1:2,5,8,10,13,16,19)]],2,function(y) sum(!is.na(y))/dim(x)[1])))
  return(profiled)
}

orthogonal_profile = rbind(ldply(shuffled_pos,function(y) profile_set(y)),profile_set(restricted),profile_set(TE_pos))
test = profile_set(list(TE_neg))
test$`.id` = NA
orthogonal_profile = rbind(orthogonal_profile,test)
rm(test)

colnames(orthogonal_profile)[1:3] = c("State","Number","Length")
orthogonal_profile$SVA = orthogonal_profile$Other
orthogonal_profile$Other = orthogonal_profile$`DNA?`+orthogonal_profile$`LINE?` + orthogonal_profile$`LTR?` + orthogonal_profile$RC + orthogonal_profile$`SINE?` + orthogonal_profile$Unknown + orthogonal_profile$`Unknown?`
orthogonal_profile$Set = c(rep("Shuffled",50),rep("Restricted",5),rep("Constitutive",5),"Inactive")

ggplot(orthogonal_profile,aes(x=Set,y=Length,fill=State)) + geom_boxplot() + scale_fill_manual(values=all_state_colors) + ylim(0,350) + ylab("Median length (bp)") +
  scale_x_discrete(limits=c("Constitutive","Shuffled","Restricted","Inactive"),labels=c("Constitutively active","Shuffled constitutively active","Restricted activity","Inactive"))

ggplot(melt(orthogonal_profile[,c("Set","State",cohorts[c(1:2,5,8,10,13,16,19)])],id.vars=c("Set","State")),aes(x=Set,y=value,fill=State)) + geom_boxplot() + scale_fill_manual(values=all_state_colors,guide=FALSE) + ylim(0,1) +
  scale_x_discrete(limits=c("Constitutive","Shuffled","Restricted","Inactive"),labels=c("Constitutively active","Shuffled constitutively active","Restricted activity","Inactive")) + facet_grid(State~variable) + 
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + ylab("Proportion overlapping feature")

ggplot(melt(orthogonal_profile[,c("Set","State","DNA","LINE","LTR","SINE","SVA","Other")],id.vars=c("Set","State")),aes(x=Set,y=value,fill=State)) + geom_boxplot() + scale_fill_manual(values=all_state_colors,guide=FALSE) + ylim(0,1) +
  scale_x_discrete(limits=c("Constitutive","Shuffled","Restricted","Inactive"),labels=c("Constitutively active","Shuffled constitutively active","Restricted activity","Inactive")) + facet_grid(State~variable) + 
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + ylab("Proportion in class")
