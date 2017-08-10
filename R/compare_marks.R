# Correlation between chromHMM, WGBS, DNase, and H3K27ac
# See 5/16/2017, 6/7/2017, 7/31/2017

# Number of TEs x sample in each WGBS x Dnase x H3K27ac state
compare_marks_unique = read.table(file="compare_marks/combine_marks_counts.txt",sep='\t')
colnames(compare_marks_unique) = c("WGBS","DNase","H3K27ac","TE_sample")
compare_marks_unique$DNase = factor(compare_marks_unique$DNase,levels=c("yes","no"))
compare_marks_unique$H3K27ac = factor(compare_marks_unique$H3K27ac,levels=c("yes","no"))
compare_marks_unique$WGBS = factor(compare_marks_unique$WGBS,levels=c(meth_states,"Hyper","Hypo"))
compare_marks_unique[which(compare_marks_unique$WGBS == "Hyper"),]$WGBS = "Hypermethylated"
compare_marks_unique[which(compare_marks_unique$WGBS == "Hypo"),]$WGBS = "Hypomethylated"
compare_marks_unique$WGBS = factor(compare_marks_unique$WGBS,levels=meth_states)

# Number of TEs x sample in each combination
aggregate(data=compare_marks_unique,TE_sample~DNase,function(x) sum(x)/(8118555+226523729))
aggregate(data=compare_marks_unique,TE_sample~H3K27ac,function(x) sum(x)/(10934893+423092851))
aggregate(data=compare_marks_unique,TE_sample~H3K27ac+DNase,function(x) sum(x)/(1916364+4791982+2835719+185221127))
chisq.test(matrix(aggregate(data=compare_marks_unique,TE_sample~H3K27ac+DNase,sum)$TE_sample,nrow=2))
aggregate(data=compare_marks_unique,TE_sample~WGBS,function(x) sum(x)/(2219209+11847153+100619288+3730186))
aggregate(data=compare_marks_unique,TE_sample~WGBS+DNase+H3K27ac,function(x) sum(x)/54407276)

# Number of TEs x sample in each WGBS x Dnase x H3K27ac x chromHMM state
compare_marks_all = read.table(file="compare_marks/combine_marks_table_counts.txt",sep='\t')
colnames(compare_marks_all) = c("chromHMM","WGBS","DNase","H3K27ac","TE_sample")
compare_marks_all$DNase = factor(compare_marks_all$DNase,levels=c("yes","no"))
compare_marks_all$H3K27ac = factor(compare_marks_all$H3K27ac,levels=c("yes","no"))
compare_marks_all$chromHMM = factor(compare_marks_all$chromHMM,levels=chromHMM_states)
compare_marks_all$WGBS = factor(compare_marks_all$WGBS,levels=c(meth_states,"Hyper","Hypo"))
compare_marks_all[which(compare_marks_all$WGBS == "Hyper"),]$WGBS = "Hypermethylated"
compare_marks_all[which(compare_marks_all$WGBS == "Hypo"),]$WGBS = "Hypomethylated"
compare_marks_all$WGBS = factor(compare_marks_all$WGBS,levels=meth_states)
compare_marks_all = merge(compare_marks_all,compare_marks_unique,by=c("WGBS","DNase","H3K27ac"),all=TRUE)
colnames(compare_marks_all)[5:6] = c("TE_sample","TE_sample_unique")
compare_marks_all$TE_sample_percent = compare_marks_all$TE_sample/compare_marks_all$TE_sample_unique

# Number of TEs x sample in each combination
aggregate(data=compare_marks_all,TE_sample~chromHMM,function(x) sum(x)/sum(compare_marks_unique$TE_sample))
compare_marks_all[which(compare_marks_all$WGBS == "Hypomethylated" & compare_marks_all$DNase == "yes" & compare_marks_all$H3K27ac == "yes"),]