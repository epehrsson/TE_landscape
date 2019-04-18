# Select replicates
replicates = as.data.frame(matrix(c("E019","E020",
                                    "E008","E014",
                                    "E003","E016",
                                    "E055","E056",
                                    "E059","E061",
                                    "E101","E102"),ncol=2,byrow = TRUE))
colnames(replicates) = c("Sample 1","Sample 2")
replicates$Set = rep("Replicate",6)

test = as.data.frame(cbind(rep(as.vector(replicates$`Sample 1`),each=4),
                           rep(c("E032","E049","E066","E071"),times=6),
                           rep("Tissue",24)))
colnames(test) = c("Sample 1","Sample 2","Set")

replicates = rbind(replicates,test)
rm(test)

print("Add metadata to replicates")
replicates$Pair = paste(replicates$`Sample 1`,replicates$`Sample 2`,sep="/")
replicates = merge(replicates,metadata[,c("Sample","Group")],by.x="Sample 1",by.y="Sample")
replicates = merge(replicates,metadata[,c("Sample","Group")],by.x="Sample 2",by.y="Sample")
colnames(replicates)[5:6] = c("Group 1","Group 2")

# Apply over replicates, number of samples per TE x state
print("Apply over replicates")
chromHMM_replicates = apply(replicates,1,function(z) ddply(chromHMM_active_matrix[which(chromHMM_active_matrix$Sample %in% as.vector(z[1:2])),],
                                                                 .(State,chromosome,start,stop,subfamily,family,class,strand,Total),summarise,Samples=length(Sample)))
names(chromHMM_replicates) = replicates$Pair
chromHMM_replicates = ldply(chromHMM_replicates)
colnames(chromHMM_replicates)[1] = "Pair"

## Add in Group information for the pair
chromHMM_replicates = merge(chromHMM_replicates,replicates[,c("Pair","Group 1","Group 2","Set")],by="Pair")

# Add Group information to replicate matrix
chromHMM_replicates = merge(chromHMM_replicates,chromHMM_active_group,
                            by.x=c(TE_coordinates,"State","Group 1"),by.y=c(TE_coordinates,"State","Group"),all.x=TRUE)
colnames(chromHMM_replicates)[15] = "Group.Count 1"
chromHMM_replicates = merge(chromHMM_replicates,chromHMM_active_group,
                            by.x=c(TE_coordinates,"State","Group 2"),by.y=c(TE_coordinates,"State","Group"),all.x=TRUE)
colnames(chromHMM_replicates)[16] = "Group.Count 2"
chromHMM_replicates[is.na(chromHMM_replicates)] = 0
chromHMM_replicates$Group.Count = ifelse(chromHMM_replicates$`Group 1` == chromHMM_replicates$`Group 2`,chromHMM_replicates$`Group.Count 1`,
                                         chromHMM_replicates$`Group.Count 1`+chromHMM_replicates$`Group.Count 2`)

# Fraction of TEs in either sample found in both
print("Apply thresholds")
## Removing TEs that are found outside the pair (reviewer suggestion)
chromHMM_replicates_pair = ddply(chromHMM_replicates[which(chromHMM_replicates$Total == chromHMM_replicates$Samples),],
                              .(State,Pair,Samples,Set),summarise,Count=length(Samples))
chromHMM_replicates_pair$Threshold = rep("Exclusive",dim(chromHMM_replicates_pair)[1])

## Removing TEs that are in 5+ samples
chromHMM_replicates_5 = ddply(chromHMM_replicates[which(chromHMM_replicates$Total < 5),],
                              .(State,Pair,Samples,Set),summarise,Count=length(Samples))
chromHMM_replicates_5$Threshold = rep("<5 samples",dim(chromHMM_replicates_5)[1])

## Removing TEs found outside the Group(s)
chromHMM_replicates_group = ddply(chromHMM_replicates[which(chromHMM_replicates$Total == chromHMM_replicates$Group.Count),],
                                  .(State,Pair,Samples,Set),summarise,Count=length(Samples))
chromHMM_replicates_group$Threshold = rep("Group exclusive",dim(chromHMM_replicates_group)[1])

## All TEs 
chromHMM_replicates_all = ddply(chromHMM_replicates,.(State,Pair,Samples,Set),summarise,Count=length(Samples))
chromHMM_replicates_all$Threshold = rep("All",dim(chromHMM_replicates_all)[1])

## Combine all
print("Combine")
chromHMM_replicates_count = rbind(chromHMM_replicates_pair,chromHMM_replicates_5,chromHMM_replicates_group,chromHMM_replicates_all)
chromHMM_replicates_count = dcast(chromHMM_replicates_count,Pair+Set+State+Threshold~Samples,value.var="Count")
chromHMM_replicates_count[is.na(chromHMM_replicates_count)] = 0
chromHMM_replicates_count$Fraction = chromHMM_replicates_count$`2`/(chromHMM_replicates_count$`1`+chromHMM_replicates_count$`2`)
chromHMM_replicates_count$Threshold = factor(chromHMM_replicates_count$Threshold,levels=c("Exclusive","<5 samples","Group exclusive","All"))
rm(list=c("chromHMM_replicates_pair","chromHMM_replicates_5","chromHMM_replicates_group","chromHMM_replicates"))