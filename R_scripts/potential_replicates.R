# Identify replicates
## Load matrices
print("Load correlation matrices")

load("raw_data/correlation_matrices/cor_H3K27me3_13.RData")
h3k27me3 = markcor

load("raw_data/correlation_matrices/cor_H3K36me3_4.RData")
h3k36me3 = markcor

load("raw_data/correlation_matrices/cor_H3K4me1_7.RData")
h3k4me1 = markcor

load("raw_data/correlation_matrices/cor_H3K4me3_1.RData")
h3k4me3 = markcor

load("raw_data/correlation_matrices/cor_H3K9me3_9.RData")
h3k9me3 = markcor

print("Rank correlation")
marks = ls(pattern="h3")
mark_cor = lapply(marks,function(x) transform_matrix(get(x)))
names(mark_cor) = marks
mark_cor = ldply(mark_cor)
colnames(mark_cor)[1] = "Mark"

mark_cor_ranked = ddply(mark_cor,.(`Sample 1`,`Sample 2`),summarise,Rank=mean(Rank))
mark_cor_ranked = mark_cor_ranked[order(mark_cor_ranked$Rank),]

replicates = mark_cor_ranked[c(1:2,5:7,13:15,17,19,21,23:24,31,33,38,40,46:47,59),]
rm(list=c("mark_cor_ranked","mark_cor","markcor","h3k27me3","h3k36me3","h3k4me1","h3k4me3","h3k9me3","marks","transform_matrix"))

# Replicate analysis
print("Add metadata to replicates")
replicates$Pair = paste(replicates$`Sample 1`,replicates$`Sample 2`,sep="/")
replicates = merge(replicates,metadata[,c("Sample","Group")],by.x="Sample 1",by.y="Sample")
replicates = merge(replicates,metadata[,c("Sample","Group")],by.x="Sample 2",by.y="Sample")
colnames(replicates)[5:6] = c("Group 1","Group 2")
replicates = replicates[order(replicates$Rank),]

## Load TE x sample x state
print("Load state matrices")
chromHMM_active_matrix = lapply(states[c(1:3,6:7)],function(x) load_state(x))
names(chromHMM_active_matrix) = states[c(1:3,6:7)]
chromHMM_active_matrix = ldply(chromHMM_active_matrix)
colnames(chromHMM_active_matrix)[1] = "State"

## Apply over replicates, number of samples per TE x state
print("Apply over replicates")
chromHMM_replicates = apply(replicates,1,function(z) ddply(chromHMM_active_matrix[which(chromHMM_active_matrix$Sample %in% as.vector(z[1:2])),],
                                                                 .(State,chromosome,start,stop,subfamily,family,class,strand,Total),summarise,Samples=length(Sample)))
rm(chromHMM_active_matrix)
names(chromHMM_replicates) = replicates$Pair
chromHMM_replicates = ldply(chromHMM_replicates)
colnames(chromHMM_replicates)[1] = "Pair"

## Add in Group information for the pair (Sample 1 only)
chromHMM_replicates = merge(chromHMM_replicates,replicates[,c("Pair","Group 1","Group 2")],by="Pair")

### Generate matrices of TE x samples in category by state
print("Load group matrices")
#load_state_group("1_TssA")
#load_state_group("2_TssAFlnk")
#load_state_group("3_TxFlnk")
#load_state_group("6_EnhG")
#load_state_group("7_Enh")

### Combine matrices
chromHMM_active = lapply(states[c(1:3,6:7)],function(x) mget(load(paste("chromHMM/chromHMM_",x,".RData",sep="")))$state_group)
names(chromHMM_active) = states[c(1:3,6:7)]
chromHMM_active = ldply(chromHMM_active)
colnames(chromHMM_active)[1] = "State"
colnames(chromHMM_active)[11] = "Group.Count"

### Add Group information to replicate matrix
chromHMM_replicates = merge(chromHMM_replicates,
                            chromHMM_active[which(chromHMM_active$Category == "Group"),c(TE_coordinates,"State","Grouping","Group.Count")],
                            by.x=c(TE_coordinates,"State","Group 1"),by.y=c(TE_coordinates,"State","Grouping"),all.x=TRUE)
colnames(chromHMM_replicates)[14] = "Group.Count 1"
chromHMM_replicates = merge(chromHMM_replicates,
                            chromHMM_active[which(chromHMM_active$Category == "Group"),c(TE_coordinates,"State","Grouping","Group.Count")],
                            by.x=c(TE_coordinates,"State","Group 2"),by.y=c(TE_coordinates,"State","Grouping"),all.x=TRUE)
colnames(chromHMM_replicates)[15] = "Group.Count 2"
chromHMM_replicates[is.na(chromHMM_replicates)] = 0
chromHMM_replicates$Group.Count = ifelse(chromHMM_replicates$`Group 1` == chromHMM_replicates$`Group 2`,chromHMM_replicates$`Group.Count 1`,
                                         chromHMM_replicates$`Group.Count 1`+chromHMM_replicates$`Group.Count 2`)
rm(chromHMM_active)

## Fraction of TEs in either sample found in both
print("Apply thresholds")
### Removing TEs that are found outside the pair (reviewer suggestion)
chromHMM_replicates_pair = ddply(chromHMM_replicates[which(chromHMM_replicates$Total == chromHMM_replicates$Samples),],
                              .(State,Pair,Samples),summarise,Count=length(Samples))
chromHMM_replicates_pair$Threshold = rep("Exclusive",dim(chromHMM_replicates_pair)[1])

### Removing TEs that are in 5+ samples
chromHMM_replicates_5 = ddply(chromHMM_replicates[which(chromHMM_replicates$Total < 5),],
                              .(State,Pair,Samples),summarise,Count=length(Samples))
chromHMM_replicates_5$Threshold = rep("<5 samples",dim(chromHMM_replicates_5)[1])

### Removing TEs found outside the Group(s)
chromHMM_replicates_group = ddply(chromHMM_replicates[which(chromHMM_replicates$Total == chromHMM_replicates$Group.Count),],
                                  .(State,Pair,Samples),summarise,Count=length(Samples))
chromHMM_replicates_group$Threshold = rep("Group exclusive",dim(chromHMM_replicates_group)[1])

### All TEs 
chromHMM_replicates_all = ddply(chromHMM_replicates,.(State,Pair,Samples),summarise,Count=length(Samples))
chromHMM_replicates_all$Threshold = rep("All",dim(chromHMM_replicates_all)[1])

## Combine all
print("Combine")
chromHMM_replicates_count = rbind(chromHMM_replicates_pair,chromHMM_replicates_5,chromHMM_replicates_group,chromHMM_replicates_all)
chromHMM_replicates_count = dcast(chromHMM_replicates_count,Pair+State+Threshold~Samples,value.var="Count")
chromHMM_replicates_count[is.na(chromHMM_replicates_count)] = 0
chromHMM_replicates_count$Fraction = chromHMM_replicates_count$`2`/(chromHMM_replicates_count$`1`+chromHMM_replicates_count$`2`)
chromHMM_replicates_count$Threshold = factor(chromHMM_replicates_count$Threshold,levels=c("Exclusive","<5 samples","Group exclusive","All"))

# Analysis 2
## For TEs in 2+ samples, how often is it the same category?
#active_reg_multi = chromHMM_active[which(chromHMM_active$Total > 1),]
#active_reg_multi = ddply(active_reg_multi,.(Category,Grouping,Category.Total,State),summarise,
#                         Exclusive=sum(Count == Total),All=length(chromosome))
#active_reg_multi$Fraction = active_reg_multi$Exclusive/(active_reg_multi$Exclusive + active_reg_multi$Shared)

## Plot
#ggplot(active_reg_multi,aes(x=Category.Total,y=Fraction,color=State)) + geom_point() + scale_color_manual(values=all_state_colors) + 
#  xlab("Total samples in category") + ylab("Fraction of exclusive TEs / shared TEs")