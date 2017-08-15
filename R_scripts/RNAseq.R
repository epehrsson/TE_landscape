# RNA analysis (TEs and exons)
# See 3/1/2017, 3/2/2017, 3/3/2017, 3/7/2017, 6/8/2017, 6/11/2017, 6/14/2017, 8/1/2017

# Normalization factors per sample
RNA_metadata = read.table("RNA_average/all.EGID.N.readlength",sep='\t',header=TRUE)

# RNA-seq sample metadata
RNA_samples = as.data.frame(c(as.vector(read.table("RNA_average/RNA_samples_stranded.txt")$V1),as.vector(read.table("RNA_average/RNA_samples_agnostic.txt")$V1)))
colnames(RNA_samples)[1] = "File"
RNA_samples$Sample = apply(RNA_samples,1,function(x) unlist(strsplit(x,"[.]"))[1])
RNA_samples$Strand = apply(RNA_samples,1,function(x) unlist(strsplit(x,"[.]"))[2])
RNA_samples[which(RNA_samples$Strand == "wig"),]$Strand = "non"
RNA_samples$Replicate = rep(NA,116)
RNA_samples[3:22,]$Replicate = apply(RNA_samples[3:22,],1,function(x) unlist(strsplit(x[2],"_"))[2])
RNA_samples[3:22,]$Sample = apply(RNA_samples[3:22,],1,function(x) unlist(strsplit(x[2],"_"))[1])
RNA_samples = merge(RNA_samples,RNA_metadata,by.x="Sample",by.y="EGID")
RNA_samples$Factor = apply(RNA_samples,1,function(x) (1000/as.numeric(x[6]))/(as.numeric(x[5])/1000000))

# TEs
# Average expression per TE
RNA_TE = read.table("RNA_average/rmsk_TEother_average.txt",sep='\t')
colnames(RNA_TE) = c("chromosome","start","stop","subfamily","class","family","strand",as.vector(RNA_samples$File))
RNA_TE[,8:123] = t(t(RNA_TE[,8:123])*RNA_samples[match(colnames(RNA_TE[,8:123]),RNA_samples$File),]$Factor)
RNA_TE = RNA_TE[,c(1:118,121)]
RNA_TE = RNA_TE[,c(1:11,14:15,18:19,22:23,26:27,30:119)]
RNA_TE[,which(RNA_samples[match(colnames(RNA_TE),RNA_samples$File),]$Strand == "neg")] = RNA_TE[,which(RNA_samples[match(colnames(RNA_TE),RNA_samples$File),]$Strand == "neg")]*-1

# Average expression per TE, replicates removed, stranded combined
test = t(ldply(unique(RNA_samples[which(RNA_samples$File %in% colnames(RNA_TE) & RNA_samples$Strand != "non"),]$Sample),function(x) apply(RNA_TE[,as.vector(RNA_samples[which(RNA_samples$Sample == x & (RNA_samples$Replicate == "r1a" | is.na(RNA_samples$Replicate))),]$File)],1,sum)))
colnames(test) = unique(RNA_samples[which(RNA_samples$File %in% colnames(RNA_TE) & RNA_samples$Strand != "non"),]$Sample)
RNA_TE_agnostic = cbind(RNA_TE[,1:7],test,RNA_TE[,106:109])
colnames(RNA_TE_agnostic)[57:60] = c("E028","E037","E038","E062")

# Number of samples TE is expressed >1 RPKM
RNA_TE_agnostic$Expressed_samples = apply(RNA_TE_agnostic[,8:60],1,function(x) sum(x > 1))

# Max expression per TE
RNA_TE_agnostic$Max_expression = apply(RNA_TE_agnostic[,8:60],1,max)

# TE analyses
mean(apply(RNA_TE_agnostic[,8:60],1,function(x) sum(x > 0)))
mean(RNA_TE_agnostic$Expressed_samples)
mean(RNA_TE_agnostic[which(RNA_TE_agnostic$Expressed_samples > 0),]$Expressed_samples)
median(RNA_TE_agnostic$Max_expression)
median(RNA_TE_agnostic[which(RNA_TE_agnostic$Expressed_samples > 0),]$Max_expression)

# TE expression by class
RNA_TE_agnostic$class_update = RNA_TE_agnostic$class
RNA_TE_agnostic$class_update = factor(RNA_TE_agnostic$class_update,levels=c("DNA","LINE","LTR","SINE","Other","RC","Unconfident"))
RNA_TE_agnostic[which(RNA_TE_agnostic$class %in% c("DNA?","LINE?","LTR?","SINE?","Unknown","Unknown?")),]$class_update = "Unconfident"
merge(aggregate(data=RNA_TE_agnostic,Expressed_samples~class_update,function(x) sum(x > 0)/length(x)),aggregate(data=RNA_TE_agnostic,Max_expression~class_update,median),by=c("class_update"))
aggregate(data=RNA_TE_agnostic[which(RNA_TE_agnostic$Expressed_samples > 0),],Max_expression~class_update,median)

# TE expression versus age
RNA_TE_agnostic = merge(RNA_TE_agnostic,rmsk_TEother_age[,c(1:7,9)],by=c("chromosome","start","stop","strand","class","family","subfamily"))
by(RNA_TE_agnostic,RNA_TE_agnostic$class_update,function(x) unlist(cor.test(x$JC_distance,x$Expressed_samples))[c("p.value","estimate.cor")])
by(RNA_TE_agnostic,RNA_TE_agnostic$class_update,function(x) unlist(cor.test(x$JC_distance,x$Max_expression))[c("p.value","estimate.cor")])

# Average expression per subfamily
RNA_TE_agnostic_subfamily = aggregate(data=RNA_TE_agnostic[,c(4:6,8:60)],.~class+family+subfamily,mean)

# Max average expression per subfamily
RNA_TE_agnostic_subfamily$Max_expression = apply(RNA_TE_agnostic_subfamily[,4:56],1,max)

# Number of samples mean subfamily expression is >1 RPKM
RNA_TE_agnostic_subfamily$Expressed_samples = apply(RNA_TE_agnostic_subfamily[,4:56],1,function(x) sum(x > 1))
dim(RNA_TE_agnostic_subfamily[which(RNA_TE_agnostic_subfamily$Max_expression > 1),])

# TE expression by feature overlap
RNA_TE_agnostic_feature = apply(rmsk_TEother_feature[,9:16],2,function(y) apply(merge(rmsk_TEother_feature[which(y != "NA"),1:7],RNA_TE_agnostic,by=c("chromosome","start","stop","strand","class","family","subfamily"))[,61:62],2,function(x) mean(x)/0.53))
RNA_TE_agnostic_feature[2,] = RNA_TE_agnostic_feature[2,]*0.53

# Exons
# Average expression per Refseq exon
RNA_refseq_exon = read.table("RNA_average/refseq_exons_average.txt",sep='\t')
colnames(RNA_refseq_exon) = c("chromosome","start","stop","exon","X","strand",as.vector(RNA_samples$File))
RNA_refseq_exon[,7:122] = t(t(RNA_refseq_exon[,7:122])*RNA_samples[match(colnames(RNA_refseq_exon[,7:122]),RNA_samples$File),]$Factor)
RNA_refseq_exon = RNA_refseq_exon[,c(1:117,120)]
RNA_refseq_exon = RNA_refseq_exon[,c(1:10,13:14,17:18,21:22,25:26,29:118)]
RNA_refseq_exon[,which(RNA_samples[match(colnames(RNA_refseq_exon),RNA_samples$File),]$Strand == "neg")] = RNA_refseq_exon[,which(RNA_samples[match(colnames(RNA_refseq_exon),RNA_samples$File),]$Strand == "neg")]*-1

# Average expression per exon, replicates removed, stranded combined
test = t(ldply(unique(RNA_samples[which(RNA_samples$File %in% colnames(RNA_refseq_exon) & RNA_samples$Strand != "non"),]$Sample),function(x) apply(RNA_refseq_exon[,as.vector(RNA_samples[which(RNA_samples$Sample == x & (RNA_samples$Replicate == "r1a" | is.na(RNA_samples$Replicate))),]$File)],1,sum)))
colnames(test) = unique(RNA_samples[which(RNA_samples$File %in% colnames(RNA_refseq_exon) & RNA_samples$Strand != "non"),]$Sample)
RNA_refseq_exon_agnostic = cbind(RNA_refseq_exon[,1:6],test,RNA_refseq_exon[,105:108])
colnames(RNA_refseq_exon_agnostic)[56:59] = c("E028","E037","E038","E062")

# Number of samples exon is expressed >1 RPKM
RNA_refseq_exon_agnostic$Expressed_samples = apply(RNA_refseq_exon_agnostic[,7:59],1,function(x) sum(x > 1))

# Max expression per exon
RNA_refseq_exon_agnostic$Max_expression = apply(RNA_refseq_exon_agnostic[,7:59],1,max)

# Exon analyses
mean(apply(RNA_refseq_exon_agnostic[,7:59],1,function(x) sum(x > 0)))
mean(RNA_refseq_exon_agnostic$Expressed_samples)
mean(RNA_refseq_exon_agnostic[which(RNA_refseq_exon_agnostic$Expressed_samples > 0),]$Expressed_samples)
median(RNA_refseq_exon_agnostic$Max_expression)
median(RNA_refseq_exon_agnostic[which(RNA_refseq_exon_agnostic$Expressed_samples > 0),]$Max_expression)
dim(unique(RNA_refseq_exon_agnostic[,c(1:3,60:61)])) #Same when stats columns are removed
mean(apply(unique(RNA_refseq_exon_agnostic[,c(1:3,7:59)]),1,function(x) sum(as.numeric(x[4:56]) > 0)))
mean(unique(RNA_refseq_exon_agnostic[,c(1:3,60:61)])$Expressed_samples)
dim(unique(RNA_refseq_exon_agnostic[which(RNA_refseq_exon_agnostic$Expressed_samples > 0),c(1:3,60:61)]))
mean(unique(RNA_refseq_exon_agnostic[which(RNA_refseq_exon_agnostic$Expressed_samples > 0),c(1:3,60:61)])$Expressed_samples)
hist(unique(RNA_refseq_exon_agnostic[,c(1:3,60:61)])$Expressed_samples)
median(unique(RNA_refseq_exon_agnostic[,c(1:3,60:61)])$Max_expression)
median(unique(RNA_refseq_exon_agnostic[which(RNA_refseq_exon_agnostic$Expressed_samples > 0),c(1:3,60:61)])$Max_expression)
