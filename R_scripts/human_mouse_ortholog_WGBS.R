# Comparison of human-mouse TE ortholog methylation and chromHMM state
# See 5/11/2017, 5/12/2017, 5/22/2017, 5/29/2017, 5/30/2017, 7/21/2017, 7/26/2017, 7/31/2017, 8/1/2017, 

load("R_datasets/TE_meth_average.RData")
source("R_scripts/WGBS_subfamily_enrichment_TE.R")

# Corresponding human-mouse samples
human_mouse_samples = read.table("Mouse/human_mouse_samples.txt",sep='\t')
colnames(human_mouse_samples) = c("Mouse","Human")
test = rbind(expand.grid(human_mouse_samples$Mouse,human_mouse_samples$Human),expand.grid(human_mouse_samples$Mouse,human_mouse_samples$Mouse),expand.grid(human_mouse_samples$Human,human_mouse_samples$Human))
colnames(test) = colnames(human_mouse_samples)
human_mouse_samples = rbind(human_mouse_samples,test)
rm(test)
human_mouse_samples = unique(human_mouse_samples[,1:2])
human_mouse_samples$Test = c(rep("Paired",13),rep("Random_Pair",104),rep("Mouse",169),rep("Human",81))
human_mouse_samples$Pair = apply(human_mouse_samples,1,function(x) paste(x[1],x[2],sep="_"))

# mm10 TEs with WGBS average methylation
mm10_rmsk_TE_WGBS = read.table("Mouse/WGBS/mm10_rmsk_TE_WGBS_avg.txt",sep='\t',header=TRUE)

# Human TEs with mouse orthologs (mm10)
human_mouse_orthologs_mm10 = read.table("Mouse/liftover/hg19_mm10_TE_intersect_same.bed",sep='\t')
colnames(human_mouse_orthologs_mm10) = c("human_chr_mm10","human_start_mm10","human_stop_mm10","human_strand_mm10","human_chr_hg19","human_start_hg19","human_stop_hg19","human_subfamily","human_class","human_family","human_strand_hg19","mouse_chr_mm10","mouse_start_mm10","mouse_stop_mm10","mouse_subfamily","mouse_class","mouse_family","mouse_strand_mm10","overlap")

# Orthologous TEs with methylation level
human_mouse_orthologs_mm10 = merge(human_mouse_orthologs_mm10,TE_meth_average[,c(1:7,match(sort(as.vector(human_mouse_samples$Human)[1:13]),colnames(TE_meth_average)))],by.y=c("chromosome","start","stop","subfamily","class","family","strand"),by.x=c("human_chr_hg19","human_start_hg19","human_stop_hg19","human_subfamily","human_class","human_family","human_strand_hg19"),all.x=TRUE)
human_mouse_orthologs_mm10 = merge(human_mouse_orthologs_mm10,mm10_rmsk_TE_WGBS,by.y=c("chromosome","start","stop","subfamily","class","family","strand"),by.x=c("mouse_chr_mm10","mouse_start_mm10","mouse_stop_mm10","mouse_subfamily","mouse_class","mouse_family","mouse_strand_mm10"),all.x=TRUE)

# Spearman correlation between all TEs for each sample pair
human_mouse_samples$Spearman = as.numeric(apply(human_mouse_samples,1,function(x) unlist(cor.test(human_mouse_orthologs_mm10[,x[1]],human_mouse_orthologs_mm10[,x[2]],method="spearman"))["estimate.rho"]))

# Human TEs with mouse orthologs, mm10, with methylation state
human_mouse_orthologs_mm10_category = human_mouse_orthologs_mm10
human_mouse_orthologs_mm10_category[,20:32] = apply(human_mouse_orthologs_mm10[,20:32],2,function(x) ifelse(x > 0.7,"Hypermethylated",ifelse(x < 0.3,"Hypomethylated",ifelse(x <= 0.7 & x >= 0.3,"Intermediate","NA"))))
human_mouse_orthologs_mm10_category[,33:41] = apply(human_mouse_orthologs_mm10[,33:41],2,function(x) ifelse(x > 70,"Hypermethylated",ifelse(x < 30,"Hypomethylated",ifelse(x <= 70 & x >= 30,"Intermediate","NA"))))

# Subfamily analysis
# Proportion of mm10 subfamily hypomethylated
mm10_rmsk_TE_WGBS_subfamily_hypo = aggregate(data=mm10_rmsk_TE_WGBS[,c(4:6,8:16)],.~subfamily+class+family,function(x) sum(na.omit(x) < 30)/length(x),na.action=na.pass)

# Proportion of mm10 and hg19 subfamily hypomethylated
hg19_mm10_TE_WGBS_subfamily_hypo = merge(TE_meth_subfamily[["Hypomethylated"]][,c(1:3,match(sort(as.vector(human_mouse_samples[1:13,]$Human)),colnames(TE_meth_subfamily[["Hypomethylated"]])))],mm10_rmsk_TE_WGBS_subfamily_hypo[,c(1,4:12)])

# Subfamily methylation values for paired samples
hg19_mm10_TE_WGBS_subfamily_hypo_paired = apply(human_mouse_samples[1:13,],1,function(x) hg19_mm10_TE_WGBS_subfamily_hypo[,c(1:3,match(x[1],colnames(hg19_mm10_TE_WGBS_subfamily_hypo)),match(x[2],colnames(hg19_mm10_TE_WGBS_subfamily_hypo)))])
names(hg19_mm10_TE_WGBS_subfamily_hypo_paired) = human_mouse_samples[1:13,]$Pair
hg19_mm10_TE_WGBS_subfamily_hypo_paired = lapply(hg19_mm10_TE_WGBS_subfamily_hypo_paired,function(x) {colnames(x)[4:5] <- c("Mouse","Human"); x})
hg19_mm10_TE_WGBS_subfamily_hypo_paired = ldply(hg19_mm10_TE_WGBS_subfamily_hypo_paired)
colnames(hg19_mm10_TE_WGBS_subfamily_hypo_paired)[1] = "Pair"

# Spearman correlation between all TE subfamilies for each sample pair
human_mouse_samples$Subfamily_Spearman = as.numeric(apply(human_mouse_samples,1,function(x) unlist(cor.test(hg19_mm10_TE_WGBS_subfamily_hypo[,x[1]],hg19_mm10_TE_WGBS_subfamily_hypo[,x[2]],method="spearman"))["estimate.rho"]))

# Add chromHMM
# TE orthologs hypomethylated in both human and mouse
human_mouse_orthologs_mm10_hypo = apply(human_mouse_samples[1:13,],1,function(x) human_mouse_orthologs_mm10_category[which(human_mouse_orthologs_mm10_category[,x[1]] == "Hypomethylated" & human_mouse_orthologs_mm10_category[,x[2]] == "Hypomethylated"),1:14])
names(human_mouse_orthologs_mm10_hypo) = human_mouse_samples[1:13,]$Pair
human_mouse_orthologs_mm10_hypo = ldply(human_mouse_orthologs_mm10_hypo)
colnames(human_mouse_orthologs_mm10_hypo)[1] = "Pair"
human_mouse_orthologs_mm10_hypo$Human_sample = llply(human_mouse_orthologs_mm10_hypo$Pair,function(x) strsplit(x,"_")[[1]][2])
human_mouse_orthologs_mm10_hypo$Mouse_sample_mm10 = llply(human_mouse_orthologs_mm10_hypo$Pair,function(x) strsplit(x,"_")[[1]][1])
human_mouse_orthologs_mm10_hypo$Mouse_sample_mm9 = human_mouse_orthologs_mm10_hypo$Mouse_sample_mm10
human_mouse_orthologs_mm10_hypo$Mouse_sample_mm9 = llply(human_mouse_orthologs_mm10_hypo$Mouse_sample_mm9,function(x) ifelse(x %in% c("ENCFF592SYK","ENCFF796POX","ENCFF446IJS"),"ENCFF640YBG",ifelse(x %in% c("ENCFF895RRK","ENCFF644XMA"),"ENCFF715UGP",ifelse(x == "ENCFF467UEZ","ENCFF880AMI",ifelse(x == "ENCFF770AUO","ENCFF539POM","NA")))))
human_mouse_orthologs_mm10_hypo = as.data.frame(lapply(human_mouse_orthologs_mm10_hypo,unlist))

# Human chromHMM state for hg19 TEs hypomethylated in mouse/human
write.table(unique(human_mouse_orthologs_mm10_hypo[,9:16]),file="Mouse/hg19_ortholog_hypo.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
hg19_orthologs_hypo_chromHMM = read.table("Mouse/hg19_ortholog_hypo_chromHMM.txt",sep='\t')[,1:9]
colnames(hg19_orthologs_hypo_chromHMM) = c(colnames(human_mouse_orthologs_mm10_hypo)[c(9:16)],"Human_state")

# mm9 chromHMM state for mm10 TEs hypomethylated in mouse/humans
write.table(unique(human_mouse_orthologs_mm10_hypo[,2:8]),file="Mouse/mm10_ortholog_hypo.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
mm10_orthologs_hypo_chromHMM = read.table("Mouse/mm10_ortholog_hypo_chromHMM.txt",sep='\t')[,1:9]
colnames(mm10_orthologs_hypo_chromHMM) = c(colnames(human_mouse_orthologs_mm10_hypo)[2:8],"Mouse_State","Mouse_sample_mm9")

# Combine chromHMM and WGBS
human_mouse_orthologs_mm10_hypo = merge(human_mouse_orthologs_mm10_hypo,hg19_orthologs_hypo_chromHMM,by=colnames(human_mouse_orthologs_mm10_hypo)[9:16],all.x=TRUE)
human_mouse_orthologs_mm10_hypo = merge(human_mouse_orthologs_mm10_hypo,mm10_orthologs_hypo_chromHMM,by=colnames(human_mouse_orthologs_mm10_hypo)[c(10:16,18)],all.x=TRUE)
human_mouse_orthologs_mm10_hypo$Mouse_State = factor(human_mouse_orthologs_mm10_hypo$Mouse_State)
