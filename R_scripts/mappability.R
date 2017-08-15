# Mappability of each TE, TE subfamily
# See 11/4/2016, 11/7/2016, 2/6/2017, 2/9/2017, 3/8/2017, 5/16/2017

# Average mappability per TE (needs potential_TE_state, potential_other_state)
rmsk_TE_map36 = read.table("rmsk_TE_mappability_36mer.txt",sep='\t')
colnames(rmsk_TE_map36) = c("chromosome","start","stop","subfamily","mappability")
rmsk_TE_map36 = merge(rmsk_TE_map36,potential_TE_state,by=c("chromosome","start","stop","subfamily"))
rmsk_TE_map36 = rmsk_TE_map36[,c(1:4,7,6,8,5,9:23)]

rmsk_other_map36 = read.table("other/rmsk_other_mappability_36mer.txt",sep='\t')
colnames(rmsk_other_map36) = c("chromosome","start","stop","subfamily","mappability")
rmsk_other_map36 = merge(rmsk_other_map36,potential_other_state,by=c("chromosome","start","stop","subfamily"))
rmsk_other_map36 = rmsk_other_map36[,c(1:4,7,6,8,5,9:23)]

rmsk_TEother_map36 = rbind(rmsk_TE_map36,rmsk_other_map36)
rmsk_TEother_map36$class_update = rmsk_TEother_map36$class
rmsk_TEother_map36$class_update = factor(rmsk_TEother_map36$class_update,levels=c("DNA","LINE","SINE","LTR","RC","Other","Unconfident"))
rmsk_TEother_map36[which(rmsk_TEother_map36$class %in% c("DNA?","LINE?","SINE?","LTR?","Unknown","Unknown?")),]$class_update = "Unconfident"

# Average mappability by class
aggregate(data=rmsk_TE_map36,mappability~class,mean)
aggregate(data=rmsk_other_map36,mappability~class,mean)
mean(rmsk_other_map36[which(rmsk_other_map36$class %in% c("DNA?","LINE?","SINE?","LTR?","Unknown","Unknown?")),]$mappability)

# Correlation between mappability, % of samples Quiescent
cor(rmsk_TE_map36$mappability,rmsk_TE_map36$X15_Quies)
cor.test(rbind(rmsk_TE_map36,rmsk_other_map36)$mappability,rbind(rmsk_TE_map36,rmsk_other_map36)$X15_Quies)

# Mappability of TEs always Quiescent vs not
test = potential_TE_state[which(apply(potential_TE_state[,8:21],1,function(x) sum(x)) == 0),]
test = merge(test[,1:7],rmsk_TE_map36[,1:8])
t.test(test$mappability,rmsk_TE_map36$mappability)

test = potential_TEother_state[which(apply(potential_TEother_state[,8:21],1,function(x) sum(x)) == 0),]
test = merge(test[,1:7],rbind(rmsk_TE_map36,rmsk_other_map36)[,1:8])
t.test(test$mappability,rbind(rmsk_TE_map36,rmsk_other_map36)$mappability)

# Average mappability by subfamily
rmsk_TE_map36_subfamily = aggregate(data=rmsk_TE_map36,mappability~subfamily+family+class,mean)
rmsk_other_map36_subfamily = aggregate(data=rmsk_other_map36,mappability~subfamily+family+class,mean)
rmsk_TEother_map36_subfamily = rbind(rmsk_TE_map36_subfamily,rmsk_other_map36_subfamily)

rmsk_TEother_map36_subfamily$class_update = rmsk_TEother_map36_subfamily$class
rmsk_TEother_map36_subfamily$class_update = factor(rmsk_TEother_map36_subfamily$class_update,levels=c("DNA","LINE","SINE","LTR","RC","Other","Unconfident"))
rmsk_TEother_map36_subfamily[which(rmsk_TEother_map36_subfamily$class %in% c("DNA?","LINE?","SINE?","LTR?","Unknown","Unknown?")),]$class_update = "Unconfident"
