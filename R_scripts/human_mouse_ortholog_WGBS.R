# Comparison of human-mouse TE ortholog methylation state
# See 5/11/2017, 5/12/2017, 5/22/2017, 5/29/2017, 5/30/2017, 7/21/2017, 7/26/2017, 7/31/2017, 8/1/2017

#load("R_datasets/TE_meth_average.RData")
#source("R_scripts/WGBS_subfamily_enrichment_TE.R")

# Corresponding human-mouse samples
human_mouse_samples = read.table("Mouse/human_mouse_samples.txt",sep='\t')
colnames(human_mouse_samples) = c("Mouse","Human")
human_mouse_samples$Pair = apply(human_mouse_samples,1,function(x) paste(x[1],x[2],sep="_"))

# Human TEs with mouse orthologs (mm10)
human_mouse_orthologs_mm10 = read.table("Mouse/liftover/hg19_mm10_TE_intersect_same.bed",sep='\t')
colnames(human_mouse_orthologs_mm10) = c("human_chr_mm10","human_start_mm10","human_stop_mm10","human_strand_mm10",hg19_coordinates,mm10_coordinates,"overlap")
human_mouse_orthologs_mm10$human_TE_mm10 = apply(human_mouse_orthologs_mm10,1,function(x) paste(x[1],x[2],x[3],x[4],sep="_"))
human_mouse_orthologs_mm10$human_TE_hg19 = apply(human_mouse_orthologs_mm10,1,function(x) paste(x[5],x[6],x[7],x[8],x[9],x[10],x[11],sep="_"))
human_mouse_orthologs_mm10$mouse_TE_mm10 = apply(human_mouse_orthologs_mm10,1,function(x) paste(x[12],x[13],x[14],x[15],x[16],x[17],x[18],sep="_"))

# mm10 TEs with WGBS average methylation
mm10_rmsk_TE_WGBS = read.table("Mouse/WGBS/mm10_rmsk_TE_WGBS_avg.txt",sep='\t',header=TRUE)
mm10_rmsk_TE_CpG_count = read.table("Mouse/WGBS/intersect/TEs/mm10_rmsk_TE_CpG_count.txt",sep='\t')
colnames(mm10_rmsk_TE_CpG_count) = c(TE_coordinates[c(1:4,6,5,7)],"CpGs")
mm10_rmsk_TE_CpG_count$CpGs = mm10_rmsk_TE_CpG_count$CpGs/2
mm10_rmsk_TE_WGBS = merge(mm10_rmsk_TE_WGBS,mm10_rmsk_TE_CpG_count[,TE_coordinates],by=TE_coordinates)

# Orthologous TEs with methylation level (only those with CpGs in both)
human_mouse_orthologs_WGBS = merge(human_mouse_orthologs_mm10[,c("human_chr_mm10","human_start_mm10","human_stop_mm10","human_strand_mm10",hg19_coordinates,mm10_coordinates,"overlap")],
                                   TE_meth_average[,c(TE_coordinates,as.vector(human_mouse_samples$Human))],by.y=TE_coordinates[c(1:4,6,5,7)],by.x=hg19_coordinates)
human_mouse_orthologs_WGBS = merge(human_mouse_orthologs_WGBS,mm10_rmsk_TE_WGBS,by.y=TE_coordinates[c(1:4,6,5,7)],by.x=mm10_coordinates)

# Human TEs with mouse orthologs, mm10, with methylation state
human_mouse_orthologs_WGBS_category = human_mouse_orthologs_WGBS
human_mouse_orthologs_WGBS_category[,as.vector(human_mouse_samples$Human)] = apply(human_mouse_orthologs_WGBS[,as.vector(human_mouse_samples$Human)],2,function(x) 
  ifelse(x > 0.7,"Hypermethylated",ifelse(x < 0.3,"Hypomethylated",ifelse(x <= 0.7 & x >= 0.3,"Intermediate","NA"))))
human_mouse_orthologs_WGBS_category[,as.vector(unique(human_mouse_samples$Mouse))] = apply(human_mouse_orthologs_WGBS[,as.vector(unique(human_mouse_samples$Mouse))],2,function(x) 
  ifelse(x > 70,"Hypermethylated",ifelse(x < 30,"Hypomethylated",ifelse(x <= 70 & x >= 30,"Intermediate","NA"))))

# TE orthologs hypomethylated in both human and mouse
human_mouse_orthologs_mm10_hypo = apply(human_mouse_samples,1,function(x) 
  human_mouse_orthologs_WGBS_category[which(human_mouse_orthologs_WGBS_category[,x[1]] == "Hypomethylated" & human_mouse_orthologs_WGBS_category[,x[2]] == "Hypomethylated"),c(mm10_coordinates,hg19_coordinates)])
names(human_mouse_orthologs_mm10_hypo) = human_mouse_samples$Pair
human_mouse_orthologs_mm10_hypo = ldply(human_mouse_orthologs_mm10_hypo)
colnames(human_mouse_orthologs_mm10_hypo)[1] = "Pair"
human_mouse_orthologs_mm10_hypo$Human_sample = llply(human_mouse_orthologs_mm10_hypo$Pair,function(x) strsplit(x,"_")[[1]][2])
human_mouse_orthologs_mm10_hypo$Mouse_sample_mm10 = llply(human_mouse_orthologs_mm10_hypo$Pair,function(x) strsplit(x,"_")[[1]][1])
human_mouse_orthologs_mm10_hypo$Mouse_sample_mm9 = llply(human_mouse_orthologs_mm10_hypo$Mouse_sample_mm10,function(x) ifelse(x %in% c("ENCFF592SYK","ENCFF796POX","ENCFF446IJS"),"ENCFF640YBG",
                                                                                                                             ifelse(x %in% c("ENCFF895RRK","ENCFF644XMA"),"ENCFF715UGP",
                                                                                                                                    ifelse(x == "ENCFF467UEZ","ENCFF880AMI",
                                                                                                                                           ifelse(x == "ENCFF770AUO","ENCFF539POM","NA")))))
human_mouse_orthologs_mm10_hypo = as.data.frame(lapply(human_mouse_orthologs_mm10_hypo,unlist))


# Subfamily analysis
# Proportion of mm10 subfamily hypomethylated
mm10_rmsk_TE_WGBS_subfamily_hypo = aggregate(data=mm10_rmsk_TE_WGBS[,c("subfamily",as.vector(unique(human_mouse_samples$Mouse)))],.~subfamily,function(x) sum(na.omit(x) < 30)/length(x),na.action=na.pass)

# Proportion of mm10 and hg19 subfamily hypomethylated
hg19_mm10_TE_WGBS_subfamily_hypo = merge(TE_meth_subfamily[["Hypomethylated"]][,c("subfamily","family","class_update",as.vector(human_mouse_samples$Human))],
                                         mm10_rmsk_TE_WGBS_subfamily_hypo[,c("subfamily",as.vector(human_mouse_samples$Human))])

# Subfamily methylation values for paired samples
hg19_mm10_TE_WGBS_subfamily_hypo_paired = apply(human_mouse_samples,1,function(x) hg19_mm10_TE_WGBS_subfamily_hypo[,c("subfamily","family","class_update",x[1],x[2])])
names(hg19_mm10_TE_WGBS_subfamily_hypo_paired) = human_mouse_samples$Pair
hg19_mm10_TE_WGBS_subfamily_hypo_paired = lapply(hg19_mm10_TE_WGBS_subfamily_hypo_paired,function(x) {colnames(x)[4:5] <- c("Mouse","Human"); x})
hg19_mm10_TE_WGBS_subfamily_hypo_paired = ldply(hg19_mm10_TE_WGBS_subfamily_hypo_paired)
colnames(hg19_mm10_TE_WGBS_subfamily_hypo_paired)[1] = "Pair"

# Expanded human-mouse samples
a = expand.grid(human_mouse_samples$Mouse,human_mouse_samples$Human)
b = expand.grid(human_mouse_samples$Mouse,human_mouse_samples$Mouse)
c = expand.grid(human_mouse_samples$Human,human_mouse_samples$Human)
test = rbind(a,b,c)
colnames(test) = c("Mouse","Human")
human_mouse_samples_expand = rbind(human_mouse_samples[,c("Human","Mouse")],test)
human_mouse_samples_expand$Test = c(rep("Paired",dim(human_mouse_samples)[1]),
                                    rep("Random_Pair",dim(a)[1]),
                                    rep("Mouse",dim(b)[1]),
                                    rep("Human",dim(c)[1]))
human_mouse_samples_expand = unique(human_mouse_samples_expand)
human_mouse_samples_expand$Pair = apply(human_mouse_samples_expand,1,function(x) paste(x[1],x[2],sep="_"))
rm(list=c("test","a","b","c"))

# Spearman correlation between all TEs for each sample pair
human_mouse_samples_expand$Spearman = as.numeric(apply(human_mouse_samples_expand,1,function(x) 
  unlist(cor.test(as.numeric(human_mouse_orthologs_WGBS[,x[1]]),as.numeric(human_mouse_orthologs_WGBS[,x[2]]),method="spearman"))["estimate.rho"]))

# Spearman correlation between all TE subfamilies for each sample pair
human_mouse_samples_expand$Subfamily_Spearman = as.numeric(apply(human_mouse_samples_expand,1,function(x) 
  unlist(cor.test(as.numeric(hg19_mm10_TE_WGBS_subfamily_hypo[,x[1]]),as.numeric(hg19_mm10_TE_WGBS_subfamily_hypo[,x[2]]),method="spearman"))["estimate.rho"]))
