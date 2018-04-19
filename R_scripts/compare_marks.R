# Correlation between chromHMM, WGBS, DNase, and H3K27ac
# See 1/20/2017, 2/9/2017, 5/16/2017, 6/7/2017, 7/31/2017

# Number of TEs x sample in each WGBS x Dnase x H3K27ac state
compare_marks_unique = read.table(file="compare_marks/combine_marks_counts_RNA.txt",sep='\t')
colnames(compare_marks_unique) = c("WGBS","DNase","H3K27ac","RNA","TE_sample")
compare_marks_unique$DNase = factor(compare_marks_unique$DNase,levels=c("yes","no"))
compare_marks_unique$H3K27ac = factor(compare_marks_unique$H3K27ac,levels=c("yes","no"))
compare_marks_unique$RNA = factor(compare_marks_unique$RNA,levels=c("yes","no"))
compare_marks_unique$WGBS = factor(compare_marks_unique$WGBS,levels=c(meth_states,"Hyper","Hypo"))
compare_marks_unique[which(compare_marks_unique$WGBS == "Hyper"),]$WGBS = "Hypermethylated"
compare_marks_unique[which(compare_marks_unique$WGBS == "Hypo"),]$WGBS = "Hypomethylated"
compare_marks_unique$WGBS = factor(compare_marks_unique$WGBS,levels=meth_states)

# Number of TEs x sample in each WGBS x Dnase x H3K27ac x chromHMM state
compare_marks_all = read.table(file="compare_marks/combine_marks_table_counts_RNA.txt",sep='\t')
colnames(compare_marks_all) = c("chromHMM","WGBS","DNase","H3K27ac","RNA","TE_sample")
compare_marks_all$DNase = factor(compare_marks_all$DNase,levels=c("yes","no"))
compare_marks_all$H3K27ac = factor(compare_marks_all$H3K27ac,levels=c("yes","no"))
compare_marks_all$RNA = factor(compare_marks_all$RNA,levels=c("yes","no"))
compare_marks_all$chromHMM = factor(compare_marks_all$chromHMM,levels=chromHMM_states)
compare_marks_all$WGBS = factor(compare_marks_all$WGBS,levels=c(meth_states,"Hyper","Hypo"))
compare_marks_all[which(compare_marks_all$WGBS == "Hyper"),]$WGBS = "Hypermethylated"
compare_marks_all[which(compare_marks_all$WGBS == "Hypo"),]$WGBS = "Hypomethylated"
compare_marks_all$WGBS = factor(compare_marks_all$WGBS,levels=meth_states)