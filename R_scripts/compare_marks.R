# Correlation between chromHMM, WGBS, DHS, H3K27ac, and RNA-seq annotations over individual TEs
# Includes chromHMM states (can be more than one per TE), methylation states (average methylation), 
# overlap with DHS or H3K27ac peak summits (true/false), and RPKM > 1 (true/false).
# If a sample does not have data for that technique, it is listed as "NA"

## compare_marks_unique - number of unique TE x sample instances with each combination of states, excluding chromHMM
## compare_marks_all - number of TE x sample instances with each combination of states, including chromHMM (multiple states per TE)
## WGBS_norm, DNase_norm, H3K27ac_norm, RNA_norm - number of unique TE x sample instances annotated with each state, by mark, excluding chromHMM
## by_chromHMM, by_WGBS, by_DNase, by_H3K27ac, by_RNA - tables by epigenetic mark, listing the number/proportion of TEs annotated with one state
## that are also annotated with another state generated from a different epigenetic mark

# Loads a table of the number of TE x sample instances with each combination of states, excluding chromHMM
# Because each TE can only have one annotation per mark, sums to the total number of TE x sample instances
compare_marks_unique = read.table(file="compare_marks/true_summit/combine_marks_counts_unique.txt",sep='\t',header=TRUE,na.strings = "missing")[,2:6]
colnames(compare_marks_unique) = c("H3K27ac","DNase","RNA","WGBS","TE_sample")
compare_marks_unique$DNase = factor(compare_marks_unique$DNase,levels=c("True","False"))
compare_marks_unique$H3K27ac = factor(compare_marks_unique$H3K27ac,levels=c("True","False"))
compare_marks_unique$RNA = factor(compare_marks_unique$RNA,levels=c("True","False"))
compare_marks_unique$WGBS = factor(compare_marks_unique$WGBS,levels=meth_states)

## Adds a row for TE x sample instances where the sample has only chromHMM data (all samples have chrY)
compare_marks_unique[nrow(compare_marks_unique)+1,] = list(NA,NA,NA,NA,NUM_TE*
                                                             dim(metadata[which(is.na(metadata$WGBS) & is.na(metadata$DNase) & is.na(metadata$H3K27ac) & is.na(metadata$RNA)),])[1])

## Removes chrY TEs that have RNA expression in samples without chrY (currently hardcoded)
compare_marks_unique[which(is.na(compare_marks_unique$H3K27ac) & is.na(compare_marks_unique$DNase) & is.na(compare_marks_unique$WGBS) & compare_marks_unique$RNA == "True"),]$TE_sample =
  compare_marks_unique[which(is.na(compare_marks_unique$H3K27ac) & is.na(compare_marks_unique$DNase) & is.na(compare_marks_unique$WGBS) & compare_marks_unique$RNA == "True"),]$TE_sample - 14
compare_marks_unique[which(is.na(compare_marks_unique$H3K27ac) & is.na(compare_marks_unique$DNase) & is.na(compare_marks_unique$WGBS) & compare_marks_unique$RNA == "False"),]$TE_sample =
  compare_marks_unique[which(is.na(compare_marks_unique$H3K27ac) & is.na(compare_marks_unique$DNase) & is.na(compare_marks_unique$WGBS) & compare_marks_unique$RNA == "False"),]$TE_sample - 126306

# Loads a table of the number of TE x sample instances with each combination of states, including chromHMM
# TEs can be counted multiple times if they are annotated with multiple chromHMM states in a single sample
compare_marks_all = read.table(file="compare_marks/true_summit/combine_marks_counts.txt",sep='\t',header=TRUE,na.strings = "missing")[,2:7]
colnames(compare_marks_all) = c("chromHMM","H3K27ac","DNase","RNA","WGBS","TE_sample")
compare_marks_all$DNase = factor(compare_marks_all$DNase,levels=c("True","False"))
compare_marks_all$H3K27ac = factor(compare_marks_all$H3K27ac,levels=c("True","False"))
compare_marks_all$RNA = factor(compare_marks_all$RNA,levels=c("True","False"))
compare_marks_all$WGBS = factor(compare_marks_all$WGBS,levels=meth_states)
compare_marks_all$chromHMM = factor(compare_marks_all$chromHMM,levels=chromHMM_states)

## Removes chrY TEs that have RNA expression in samples without chrY (currently hardcoded)
compare_marks_all = compare_marks_all[1:(nrow(compare_marks_all)-2),]

# Counts the number of unique TE x sample instances annotated with each state (excluding chromHMM)
WGBS_norm = aggregate(data=compare_marks_unique,TE_sample~WGBS,function(x) sum(x))
DNase_norm = aggregate(data=compare_marks_unique,TE_sample~DNase,function(x) sum(x))
H3K27ac_norm = aggregate(data=compare_marks_unique,TE_sample~H3K27ac,function(x) sum(x))
RNA_norm = aggregate(data=compare_marks_unique,TE_sample~RNA,function(x) sum(x))

# Creates tables by epigenetic mark, listing the number/proportion of TEs annotated with one state
# that are also annotated with another state generated from a different epigenetic mark.
# For comparison of non-chromHMM marks to chromHMM, 
# normalizes by the number of unique TE x sample instances in the non-chromHMM state

## chromHMM
by_chromHMM = melt(compare_marks_all,id.vars=c("chromHMM","TE_sample"),variable.name="Mark",value.name="State")
by_chromHMM = ddply(na.omit(by_chromHMM),.(chromHMM,Mark,State),summarise,TE_sample=sum(TE_sample))
by_chromHMM = ddply(by_chromHMM,.(chromHMM,Mark),transform,Proportion=TE_sample/sum(TE_sample))
by_chromHMM$State = factor(by_chromHMM$State,levels=c(chromHMM_states,meth_states,"True","False"))
by_chromHMM$Mark = factor(by_chromHMM$Mark,level=c("chromHMM","WGBS","DNase","H3K27ac","RNA"))

## WGBS
by_WGBS = rbind(melt(compare_marks_unique[which(!is.na(compare_marks_unique$WGBS)),],id.vars=c("WGBS","TE_sample")),
                melt(compare_marks_all[which(!is.na(compare_marks_all$WGBS)),c("WGBS","chromHMM","TE_sample")],id.vars=c("WGBS","TE_sample")))    
colnames(by_WGBS)[3:4] = c("Mark","State")
by_WGBS = ddply(na.omit(by_WGBS),.(WGBS,Mark,State),summarise,TE_sample=sum(TE_sample))
by_WGBS = ddply(by_WGBS,.(WGBS,Mark),transform,Proportion=TE_sample/sum(TE_sample))
by_WGBS[which(by_WGBS$Mark == "chromHMM"),]$Proportion = by_WGBS[which(by_WGBS$Mark == "chromHMM"),]$TE_sample/WGBS_norm[match(by_WGBS[which(by_WGBS$Mark == "chromHMM"),]$WGBS,WGBS_norm$WGBS),]$TE_sample
by_WGBS$State = factor(by_WGBS$State,levels=c(chromHMM_states,meth_states,"True","False"))
by_WGBS$Mark = factor(by_WGBS$Mark,level=c("chromHMM","WGBS","DNase","H3K27ac","RNA"))

## DHS
by_DNase = rbind(melt(compare_marks_unique[which(!is.na(compare_marks_unique$DNase)),],id.vars=c("DNase","TE_sample")),
                melt(compare_marks_all[which(!is.na(compare_marks_all$DNase)),c("DNase","chromHMM","TE_sample")],id.vars=c("DNase","TE_sample")))    
colnames(by_DNase)[3:4] = c("Mark","State")
by_DNase = ddply(na.omit(by_DNase),.(DNase,Mark,State),summarise,TE_sample=sum(TE_sample))
by_DNase = ddply(by_DNase,.(DNase,Mark),transform,Proportion=TE_sample/sum(TE_sample))
by_DNase[which(by_DNase$Mark == "chromHMM"),]$Proportion = by_DNase[which(by_DNase$Mark == "chromHMM"),]$TE_sample/DNase_norm[match(by_DNase[which(by_DNase$Mark == "chromHMM"),]$DNase,DNase_norm$DNase),]$TE_sample
by_DNase$State = factor(by_DNase$State,levels=c(chromHMM_states,meth_states,"True","False"))
by_DNase$Mark = factor(by_DNase$Mark,level=c("chromHMM","WGBS","DNase","H3K27ac","RNA"))

## H3K27ac
by_H3K27ac = rbind(melt(compare_marks_unique[which(!is.na(compare_marks_unique$H3K27ac)),],id.vars=c("H3K27ac","TE_sample")),
                melt(compare_marks_all[which(!is.na(compare_marks_all$H3K27ac)),c("H3K27ac","chromHMM","TE_sample")],id.vars=c("H3K27ac","TE_sample")))    
colnames(by_H3K27ac)[3:4] = c("Mark","State")
by_H3K27ac = ddply(na.omit(by_H3K27ac),.(H3K27ac,Mark,State),summarise,TE_sample=sum(TE_sample))
by_H3K27ac = ddply(by_H3K27ac,.(H3K27ac,Mark),transform,Proportion=TE_sample/sum(TE_sample))
by_H3K27ac[which(by_H3K27ac$Mark == "chromHMM"),]$Proportion = by_H3K27ac[which(by_H3K27ac$Mark == "chromHMM"),]$TE_sample/H3K27ac_norm[match(by_H3K27ac[which(by_H3K27ac$Mark == "chromHMM"),]$H3K27ac,H3K27ac_norm$H3K27ac),]$TE_sample
by_H3K27ac$State = factor(by_H3K27ac$State,levels=c(chromHMM_states,meth_states,"True","False"))
by_H3K27ac$Mark = factor(by_H3K27ac$Mark,level=c("chromHMM","WGBS","DNase","H3K27ac","RNA"))

## RNA
by_RNA = rbind(melt(compare_marks_unique[which(!is.na(compare_marks_unique$RNA)),],id.vars=c("RNA","TE_sample")),
                melt(compare_marks_all[which(!is.na(compare_marks_all$RNA)),c("RNA","chromHMM","TE_sample")],id.vars=c("RNA","TE_sample")))    
colnames(by_RNA)[3:4] = c("Mark","State")
by_RNA = ddply(na.omit(by_RNA),.(RNA,Mark,State),summarise,TE_sample=sum(TE_sample))
by_RNA = ddply(by_RNA,.(RNA,Mark),transform,Proportion=TE_sample/sum(TE_sample))
by_RNA[which(by_RNA$Mark == "chromHMM"),]$Proportion = by_RNA[which(by_RNA$Mark == "chromHMM"),]$TE_sample/RNA_norm[match(by_RNA[which(by_RNA$Mark == "chromHMM"),]$RNA,RNA_norm$RNA),]$TE_sample
by_RNA$State = factor(by_RNA$State,levels=c(chromHMM_states,meth_states,"True","False"))
by_RNA$Mark = factor(by_RNA$Mark,level=c("chromHMM","WGBS","DNase","H3K27ac","RNA"))