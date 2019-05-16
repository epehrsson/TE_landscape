# TEs
# Matrix of state x TE x Sample - all TEs ever in the state
test = dcast(chromHMM_active_matrix,State+chromosome+start+stop+subfamily+family+class+strand~Sample,length)
ddply(test,.(State),summarise,TEs=length(State))

# For each sample pair, fraction concordant
concordance = apply(replicates,1,function(z) ddply(test[,c("State",z[1:2])],.(State),function(x) sum(x[2] == x[3])/dim(x)[1]))
names(concordance) = replicates$Pair
concordance = ldply(concordance)
colnames(concordance) = c("Pair","State","Fraction")
concordance = merge(concordance,replicates[,c("Pair","Set")],by="Pair")

by(concordance,concordance$State,function(x) wilcox.test(Fraction~Set,data=x))
