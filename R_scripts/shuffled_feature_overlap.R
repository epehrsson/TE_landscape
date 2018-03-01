# Feature overlap control with shuffled TEs

load("R_datasets/shuffled_TEs_feature.RData")

shuffled_colnames = colnames(rmsk_TE_shuffled[[1]])
rmsk_TE_shuffled = lapply(rmsk_TE_shuffled,function(x) transform(x,class_update = convert_class(x$class)))
rmsk_TE_shuffled = lapply(rmsk_TE_shuffled, setNames, nm = c(shuffled_colnames,"class_update"))

feature_proportion_shuffled = ldply(rmsk_TE_shuffled,function(y) apply(y[,cohorts[1:19]],2,function(x) length(na.omit(x))))
feature_proportion_shuffled_long = melt(as.matrix(feature_proportion_shuffled/NUM_TE))
colnames(feature_proportion_shuffled_long) = c("Iteration","Feature","Proportion")

# Class 
feature_proportion_shuffled_class = ldply(rmsk_TE_shuffled,function(z) ddply(z,~class_update,function(y) apply(y[,cohorts[1:19]],2,function(x) length(na.omit(x))/length(x))))
feature_proportion_shuffled_class$Iteration = rep(seq(1,10,1),each=6)
feature_proportion_shuffled_class_long = melt(feature_proportion_shuffled_class,id.vars=c("Iteration","class_update"))
colnames(feature_proportion_shuffled_class_long)[3:4] = c("Feature","Proportion")

# Subfamily
feature_proportion_shuffled_subfamily = ldply(rmsk_TE_shuffled,function(z) ddply(z,~subfamily,function(y) apply(y[,cohorts[1:19]],2,function(x) length(na.omit(x))/length(x))))
feature_proportion_shuffled_subfamily$Iteration = rep(seq(1,10,1),each=968)
feature_proportion_shuffled_subfamily_long = melt(feature_proportion_shuffled_subfamily,id.vars=c("Iteration","subfamily"))
colnames(feature_proportion_shuffled_subfamily_long)[3:4] = c("Feature","Proportion")