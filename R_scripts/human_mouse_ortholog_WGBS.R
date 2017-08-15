# Comparison of human-mouse TE ortholog methylation and chromHMM state
# See 5/11/2017, 5/12/2017, 5/22/2017, 5/29/2017, 5/30/2017, 7/21/2017, 7/26/2017, 7/31/2017, 8/1/2017, 

# mm10 TEs with WGBS average methylation
mm10_rmsk_TE_WGBS = read.table("Mouse/mouse_WGBS/mm10_rmsk_TE_WGBS_avg.txt",sep='\t',header=TRUE)

# Human TEs with mouse orthologs, mm10, with methylation level
human_mouse_orthologs_mm10 = read.table("hg19_mm10_TE_intersect_same.bed",sep='\t')
colnames(human_mouse_orthologs_mm10) = c("human_chr_mm10","human_start_mm10","human_stop_mm10","human_strand_mm10","human_chr_hg19","human_start_hg19","human_stop_hg19","human_subfamily","human_class","human_family","human_strand_hg19","mouse_chr_mm10","mouse_start_mm10","mouse_stop_mm10","mouse_subfamily","mouse_class","mouse_family","mouse_strand_mm10","overlap")
human_mouse_orthologs_mm10 = merge(human_mouse_orthologs_mm10,TEother_meth_wCpG_average[,c(1:7,match(sort(as.vector(human_mouse_samples$Human)),colnames(TEother_meth_wCpG_average)))],by.y=c("chromosome","start","stop","subfamily","class","family","strand"),by.x=c("human_chr_hg19","human_start_hg19","human_stop_hg19","human_subfamily","human_class","human_family","human_strand_hg19"),all.x=TRUE)
human_mouse_orthologs_mm10 = merge(human_mouse_orthologs_mm10,mm10_rmsk_TE_WGBS,by.y=c("chromosome","start","stop","subfamily","class","family","strand"),by.x=c("mouse_chr_mm10","mouse_start_mm10","mouse_stop_mm10","mouse_subfamily","mouse_class","mouse_family","mouse_strand_mm10"),all.x=TRUE)

# Updating hg19 methylation levels
human_mouse_orthologs_mm10 = human_mouse_orthologs_mm10[,c(1:19,33:41)]
human_mouse_orthologs_mm10 = merge(human_mouse_orthologs_mm10,TE_meth_average[,c(1:7,match(sort(as.vector(human_mouse_samples$Human)[1:13]),colnames(TE_meth_average)))],by.y=c("chromosome","start","stop","subfamily","class","family","strand"),by.x=c("human_chr_hg19","human_start_hg19","human_stop_hg19","human_subfamily","human_class","human_family","human_strand_hg19"),all.x=TRUE)

# Number of orthologs with CpGs
merge(TE_CpG_count,unique(human_mouse_orthologs_mm10[,1:7]),by.x=colnames(TE_CpG_count)[1:7],by.y=colnames(human_mouse_orthologs_mm10)[1:7])

# Number of CpGs per orthologs with CpGs
mean(merge(TE_CpG_count,unique(human_mouse_orthologs_mm10[,1:7]),by.x=colnames(TE_CpG_count)[1:7],by.y=colnames(human_mouse_orthologs_mm10)[1:7])$CpGs)

# Human TEs with mouse orthologs, mm10, with methylation state
human_mouse_orthologs_mm10_category = human_mouse_orthologs_mm10
human_mouse_orthologs_mm10_category[,29:41] = apply(human_mouse_orthologs_mm10[,29:41],2,function(x) ifelse(x > 0.7,"Hypermethylated",ifelse(x < 0.3,"Hypomethylated",ifelse(x <= 0.7 & x >= 0.3,"Intermediate","NA"))))
human_mouse_orthologs_mm10_category[,20:28] = apply(human_mouse_orthologs_mm10[,20:28],2,function(x) ifelse(x > 70,"Hypermethylated",ifelse(x < 30,"Hypomethylated",ifelse(x <= 70 & x >= 30,"Intermediate","NA"))))

# Chi-sq test for human-mouse pairs
apply(human_mouse_samples[1:13,],1,function(x) unlist(chisq.test(matrix(c(nrow(human_mouse_orthologs_mm10_category[which(human_mouse_orthologs_mm10_category[,x[1]] == "Hypomethylated" & human_mouse_orthologs_mm10_category[,x[2]] == "Hypomethylated"),]),nrow(human_mouse_orthologs_mm10_category[which(human_mouse_orthologs_mm10_category[,x[1]] == "Hypomethylated" & human_mouse_orthologs_mm10_category[,x[2]] %in% c("Hypermethylated","Intermediate")),]),nrow(human_mouse_orthologs_mm10_category[which(human_mouse_orthologs_mm10_category[,x[1]] %in% c("Hypermethylated","Intermediate") & human_mouse_orthologs_mm10_category[,x[2]] == "Hypomethylated"),]),nrow(human_mouse_orthologs_mm10_category[which(human_mouse_orthologs_mm10_category[,x[1]] %in% c("Hypermethylated","Intermediate") & human_mouse_orthologs_mm10_category[,x[2]] %in% c("Hypermethylated","Intermediate")),])),nrow=2)))["p.value"])

# Corresponding human-mouse samples
human_mouse_samples = read.table("Mouse/human_mouse_samples.txt",sep='\t')
colnames(human_mouse_samples) = c("Mouse","Human")
test = rbind(expand.grid(human_mouse_samples$Mouse,human_mouse_samples$Human),expand.grid(human_mouse_samples$Mouse,human_mouse_samples$Mouse),expand.grid(human_mouse_samples$Human,human_mouse_samples$Human))
colnames(test) = colnames(human_mouse_samples)
human_mouse_samples = rbind(human_mouse_samples,test)
human_mouse_samples = unique(human_mouse_samples[,1:2])
human_mouse_samples$Test = c(rep("Paired",13),rep("Random_Pair",104),rep("Mouse",169),rep("Human",81))
human_mouse_samples$Pair = apply(human_mouse_samples,1,function(x) paste(x[1],x[2],sep="_"))

# Pearson correlation between all TEs for each sample pair
human_mouse_samples$Pearson = as.numeric(apply(human_mouse_samples,1,function(x) unlist(cor.test(human_mouse_orthologs_mm10[,x[1]],human_mouse_orthologs_mm10[,x[2]]))["estimate.cor"]))

# Pearson correlation between all TE subfamilies for each sample pair
human_mouse_samples$Subfamily_Pearson = as.numeric(apply(human_mouse_samples,1,function(x) unlist(cor.test(hg19_mm10_TE_WGBS_subfamily_hypo[,x[1]],hg19_mm10_TE_WGBS_subfamily_hypo[,x[2]]))["estimate.cor"]))

# Proportion of mm10 subfamily hypomethylated
mm10_rmsk_TE_WGBS_subfamily_hypo = aggregate(data=mm10_rmsk_TE_WGBS[,c(4:6,8:16)],.~subfamily+class+family,function(x) sum(na.omit(x) < 30)/length(na.omit(x)),na.action=na.pass)

# Proportion of mm10 and hg19 subfamily hypomethylated
hg19_mm10_TE_WGBS_subfamily_hypo = merge(TE_meth_subfamily_hypo[,c(1:3,match(human_mouse_samples[1:13,]$Human,colnames(TE_meth_subfamily_hypo)))],mm10_rmsk_TE_WGBS_subfamily_hypo[,c(1,4:12)])

# Subfamily methylation values for paired samples
test = apply(human_mouse_samples[1:13,],1,function(x) hg19_mm10_TE_WGBS_subfamily_hypo[,c(1:3,match(x[1],colnames(hg19_mm10_TE_WGBS_subfamily_hypo)),match(x[2],colnames(hg19_mm10_TE_WGBS_subfamily_hypo)))])
names(test) = human_mouse_samples[1:13,]$Pair
test = lapply(test,function(x) {colnames(x)[4:5] <- c("Mouse","Human"); x})
test = ldply(test)
colnames(test)[1] = "Pair"
hg19_mm10_TE_WGBS_subfamily_hypo_paired = test

# chromHMM state for mm10 TEs hypo in humans (mm9)
human_mouse_orthologs_mm10_hypo_chromHMM = read.table("Mouse/mouse_human_ortholog_hypo_chromHMM_sum.txt",sep='\t')
colnames(human_mouse_orthologs_mm10_hypo_chromHMM) = c(colnames(human_mouse_orthologs_mm10_hypo)[10:16],"Mouse_State","Mouse_sample_mm9","Bases")

# Human/mouse TE orthologs hypomethylated in both, with chromHMM state
human_mouse_orthologs_mm10_hypo = apply(human_mouse_samples[1:13,],1,function(x) human_mouse_orthologs_mm10_category[which(human_mouse_orthologs_mm10_category[,x[1]] == "Hypomethylated" & human_mouse_orthologs_mm10_category[,x[2]] == "Hypomethylated"),1:14])
names(human_mouse_orthologs_mm10_hypo) = human_mouse_samples[1:13,]$Pair
human_mouse_orthologs_mm10_hypo = ldply(human_mouse_orthologs_mm10_hypo)
colnames(human_mouse_orthologs_mm10_hypo)[1] = "Pair"
human_mouse_orthologs_mm10_hypo$Human_sample = llply(human_mouse_orthologs_mm10_hypo$Pair,function(x) strsplit(x,"_")[[1]][2])
human_mouse_orthologs_mm10_hypo$Mouse_sample = llply(human_mouse_orthologs_mm10_hypo$Pair,function(x) strsplit(x,"_")[[1]][1])
human_mouse_orthologs_mm10_hypo$Mouse_sample_mm9 = human_mouse_orthologs_mm10_hypo$Mouse_sample
human_mouse_orthologs_mm10_hypo$Mouse_sample_mm9 = llply(human_mouse_orthologs_mm10_hypo$Mouse_sample_mm9,function(x) ifelse(x %in% c("ENCFF592SYK","ENCFF796POX","ENCFF446IJS"),"ENCFF640YBG",ifelse(x %in% c("ENCFF895RRK","ENCFF644XMA"),"ENCFF715UGP",ifelse(x == "ENCFF467UEZ","ENCFF880AMI",ifelse(x == "ENCFF770AUO","ENCFF539POM","NA")))))
human_mouse_orthologs_mm10_hypo = as.data.frame(lapply(human_mouse_orthologs_mm10_hypo,unlist))

# Mouse chromHMM state for TEs hypomethylated in mouse/human
mouse_human_orthologs_mm10_hypo_chromHMM = read.table("Mouse/mouse_human_ortholog_hypo_chromHMM_sum.txt",sep='\t')
colnames(mouse_human_orthologs_mm10_hypo_chromHMM) = c(colnames(human_mouse_orthologs_mm10_hypo)[9:15],"Mouse_State","Mouse_sample_mm9","Bases")
human_mouse_orthologs_mm10_hypo = merge(human_mouse_orthologs_mm10_hypo,mouse_human_orthologs_mm10_hypo_chromHMM[,1:9],by=colnames(human_mouse_orthologs_mm10_hypo)[c(9:15,18)],all.x=TRUE)
human_mouse_orthologs_mm10_hypo$Mouse_State = factor(human_mouse_orthologs_mm10_hypo$Mouse_State)

# Human chromHMM state for TEs hypomethylated in mouse/human
human_mouse_orthologs_mm10_chromHMM = read.table("Mouse/human_mouse_ortholog_hypo_chromHMM.txt",sep='\t')
human_mouse_orthologs_mm10_chromHMM = human_mouse_orthologs_mm10_chromHMM[,1:9]
colnames(human_mouse_orthologs_mm10_chromHMM) = c(colnames(human_mouse_orthologs_mm10_hypo)[c(10:17)],"Human_state")
human_mouse_orthologs_mm10_hypo = merge(human_mouse_orthologs_mm10_hypo,human_mouse_orthologs_mm10_chromHMM,by=colnames(human_mouse_orthologs_mm10_hypo)[10:17],all.x=TRUE)
