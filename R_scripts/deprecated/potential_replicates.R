# Apply over replicates, number of samples per TE x state
print("Apply over replicates")
chromHMM_replicates = apply(replicates,1,function(z) ddply(chromHMM_active_matrix[which(chromHMM_active_matrix$Sample %in% as.vector(z[1:2])),],
                                                                 .(State,chromosome,start,stop,subfamily,family,class,strand),summarise,Samples=length(Sample)))
names(chromHMM_replicates) = replicates$Pair
chromHMM_replicates = ldply(chromHMM_replicates)
colnames(chromHMM_replicates)[1] = "Pair"
chromHMM_replicates = merge(chromHMM_replicates,replicates[,c("Pair","Set")],by="Pair")

# Combine all active regulatory states
chromHMM_replicates_active = ddply(chromHMM_replicates,.(Pair,Set,chromosome,start,stop,subfamily,family,class,strand),summarise,Samples=max(Samples))

# Fraction of TEs in either sample found in both
chromHMM_replicates_all = ddply(chromHMM_replicates_active,.(Pair,Set,Samples),summarise,Count=length(Samples))

## Combine all
chromHMM_replicates_count = dcast(chromHMM_replicates_all,Pair+Set~Samples,value.var="Count")
chromHMM_replicates_count[is.na(chromHMM_replicates_count)] = 0
chromHMM_replicates_count$Fraction = chromHMM_replicates_count$`2`/(chromHMM_replicates_count$`1`+chromHMM_replicates_count$`2`)