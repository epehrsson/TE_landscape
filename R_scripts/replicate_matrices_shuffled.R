# Shuffled replicates

# Load TE x sample x state
## chromHMM
print("Loading active reg matrix")
chromHMM_active_shuffled = read.table("chromHMM/shuffled_TEs/rmsk_TE_shuffle_1_active_reg.txt",sep='\t')[,c(1:8,10)]
colnames(chromHMM_active_shuffled) = c(TE_coordinates[c(1:4,6,5,7)],"State","Sample")

print("Adding state total")
load("R_datasets/shuffled_chromHMM.RData")
test = shuffled_chromHMM_potential[[1]]
colnames(test)[8:22] = chromHMM_states
test = melt(test[,c(TE_coordinates,states[c(1:3,6:7)])],id.vars=TE_coordinates)
colnames(test)[8:9] = c("State","Total")
rm(shuffled_chromHMM_potential)

chromHMM_active_shuffled = merge(chromHMM_active_shuffled,test,by=c(TE_coordinates,"State"))
rm(test)

## H3K27ac
load("R_datasets/shuffled_H3K27ac.RData")
H3K27ac_pairs_shuffled = melt(shuffled_H3K27ac_potential[[1]][,c(TE_coordinates,as.vector(metadata[which(!is.na(metadata$H3K27ac)),]$Sample),"Samples")],
                              id.vars=c(TE_coordinates,"Samples"))
colnames(H3K27ac_pairs_shuffled)[8:10] = c("Total","Sample","Peaks")
rm(shuffled_H3K27ac_potential)

H3K27ac_pairs_shuffled = H3K27ac_pairs_shuffled[which(H3K27ac_pairs_shuffled$Peaks > 0),c(TE_coordinates,"Total","Sample")]
H3K27ac_pairs_shuffled$State = rep("H3K27ac",dim(H3K27ac_pairs_shuffled)[1])

chromHMM_active_shuffled = rbind(chromHMM_active_shuffled,H3K27ac_pairs_shuffled)
rm(H3K27ac_pairs_shuffled)

# Replicates
# Apply over replicates, number of samples per TE x state
print("Apply over replicates")
replicates_shuffled = apply(replicates,1,function(z) ddply(chromHMM_active_shuffled[which(chromHMM_active_shuffled$Sample %in% as.vector(z[1:2])),],
                                                           .(State,chromosome,start,stop,subfamily,family,class,strand,Total),summarise,Samples=length(Sample)))
names(replicates_shuffled) = replicates$Pair
replicates_shuffled = ldply(replicates_shuffled)
colnames(replicates_shuffled)[1] = "Pair"

## Add in Set information for the pair
replicates_shuffled = merge(replicates_shuffled,replicates[,c("Pair","Set")],by="Pair")

# Fraction of TEs in either sample found in both
print("Apply thresholds")
## Removing TEs that are found outside the pair (reviewer suggestion)
replicates_shuffled_pair = ddply(replicates_shuffled[which(replicates_shuffled$Total == replicates_shuffled$Samples),],
                                 .(State,Pair,Samples,Set),summarise,Count=length(Samples))
replicates_shuffled_pair$Threshold = rep("Exclusive",dim(replicates_shuffled_pair)[1])

## Removing TEs that are in 5+ samples
replicates_shuffled_5 = ddply(replicates_shuffled[which(replicates_shuffled$Total < 5),],
                              .(State,Pair,Samples,Set),summarise,Count=length(Samples))
replicates_shuffled_5$Threshold = rep("<5 samples",dim(replicates_shuffled_5)[1])

## All TEs 
replicates_shuffled_all = ddply(replicates_shuffled,.(State,Pair,Samples,Set),summarise,Count=length(Samples))
replicates_shuffled_all$Threshold = rep("All",dim(replicates_shuffled_all)[1])

## Combine all
print("Combine")
replicates_shuffled_count = rbind(replicates_shuffled_pair,replicates_shuffled_5,replicates_shuffled_all)
replicates_shuffled_count = dcast(replicates_shuffled_count,Pair+Set+State+Threshold~Samples,value.var="Count")
replicates_shuffled_count[is.na(replicates_shuffled_count)] = 0
replicates_shuffled_count$Fraction = replicates_shuffled_count$`2`/(replicates_shuffled_count$`1`+replicates_shuffled_count$`2`)
replicates_shuffled_count$Threshold = factor(replicates_shuffled_count$Threshold,levels=c("Exclusive","<5 samples","All"))
rm(list=c("replicates_shuffled_pair","replicates_shuffled_5","replicates_shuffled"))