# Load RNA matrices (TEs and exons)
# See 3/1/2017, 3/2/2017, 3/3/2017, 3/7/2017, 6/8/2017, 6/11/2017, 6/14/2017, 8/1/2017

# Normalization factors per sample
RNA_metadata = read.table("raw_data/RNAseq/all.EGID.N.readlength",sep='\t',header=TRUE)

# RNA-seq sample metadata
RNA_samples = as.data.frame(c(as.vector(read.table("RNAseq/RNA_samples_stranded.txt")$V1),as.vector(read.table("RNAseq/RNA_samples_agnostic.txt")$V1)))
colnames(RNA_samples)[1] = "File"
RNA_samples$Sample = apply(RNA_samples,1,function(x) unlist(strsplit(x,"[.]"))[1])
RNA_samples$Strand = apply(RNA_samples,1,function(x) unlist(strsplit(x,"[.]"))[2])
RNA_samples[which(RNA_samples$Strand == "wig"),]$Strand = "non"
RNA_samples$Replicate = rep(NA,116)
RNA_samples[3:22,]$Replicate = apply(RNA_samples[3:22,],1,function(x) unlist(strsplit(x[2],"_"))[2])
RNA_samples[3:22,]$Sample = apply(RNA_samples[3:22,],1,function(x) unlist(strsplit(x[2],"_"))[1])
RNA_samples = merge(RNA_samples,RNA_metadata,by.x="Sample",by.y="EGID")
RNA_samples$Factor = apply(RNA_samples,1,function(x) (1000/as.numeric(x[6]))/(as.numeric(x[5])/1000000))
rm(RNA_metadata)

# TEs
# Average expression per TE
RNA_TE = read.table("RNAseq/rmsk_TEother_average.txt",sep='\t')
colnames(RNA_TE) = c("chromosome","start","stop","subfamily","class","family","strand",c(as.vector(read.table("RNAseq/RNA_samples_stranded.txt")$V1),as.vector(read.table("RNAseq/RNA_samples_agnostic.txt")$V1)))
RNA_TE[,8:123] = t(t(RNA_TE[,8:123])*RNA_samples[match(colnames(RNA_TE[,8:123]),RNA_samples$File),]$Factor)
# Remove failed agnostic samples
RNA_TE = RNA_TE[,c(1:118,121)]
# Remove replicates
RNA_TE = RNA_TE[,c(1:11,14:15,18:19,22:23,26:27,30:119)]
# Combine strands
RNA_TE[,which(RNA_samples[match(colnames(RNA_TE),RNA_samples$File),]$Strand == "neg")] = RNA_TE[,which(RNA_samples[match(colnames(RNA_TE),RNA_samples$File),]$Strand == "neg")]*-1
test = t(ldply(unique(RNA_samples[which(RNA_samples$File %in% colnames(RNA_TE) & RNA_samples$Strand != "non"),]$Sample),function(x) apply(RNA_TE[,as.vector(RNA_samples[which(RNA_samples$Sample == x & (RNA_samples$Replicate == "r1a" | is.na(RNA_samples$Replicate))),]$File)],1,sum)))
colnames(test) = unique(RNA_samples[which(RNA_samples$File %in% colnames(RNA_TE) & RNA_samples$Strand != "non"),]$Sample)
RNA_TE_agnostic = cbind(RNA_TE[,1:7],test,RNA_TE[,106:109])
colnames(RNA_TE_agnostic)[57:60] = c("E028","E037","E038","E062")
rm(test)

# Number of samples TE is expressed >1 RPKM 
RNA_TE_agnostic$Expressed_samples = apply(RNA_TE_agnostic[,9:60],1,function(x) sum(x > 1))

# Number of samples TE is expressed >1 RPKM, no cancer cell lines/IMR90
RNA_TE_agnostic$Expressed_samples_noCancer = apply(RNA_TE_agnostic[,9:60],1,function(x) sum(as.numeric(x[which(metadata[match(colnames(RNA_TE_agnostic)[9:60],metadata$Sample),]$Exclude == "Include")] > 1)))

# Max expression per TE
RNA_TE_agnostic$Max_expression = apply(RNA_TE_agnostic[,9:60],1,max)

# Updating class
RNA_TE_agnostic$class_update = convert_class(RNA_TE_agnostic$class)

# Exons
# Average expression per Refseq exon
RNA_refseq_exon = read.table("RNAseq/refseq_exons_average.txt",sep='\t')
colnames(RNA_refseq_exon) = c("chromosome","start","stop","exon","X","strand",c(as.vector(read.table("RNAseq/RNA_samples_stranded.txt")$V1),as.vector(read.table("RNAseq/RNA_samples_agnostic.txt")$V1)))
RNA_refseq_exon[,7:122] = t(t(RNA_refseq_exon[,7:122])*RNA_samples[match(colnames(RNA_refseq_exon[,7:122]),RNA_samples$File),]$Factor)
RNA_refseq_exon = RNA_refseq_exon[,c(1:117,120)]
RNA_refseq_exon = RNA_refseq_exon[,c(1:10,13:14,17:18,21:22,25:26,29:118)]
RNA_refseq_exon[,which(RNA_samples[match(colnames(RNA_refseq_exon),RNA_samples$File),]$Strand == "neg")] = RNA_refseq_exon[,which(RNA_samples[match(colnames(RNA_refseq_exon),RNA_samples$File),]$Strand == "neg")]*-1

# Average expression per exon, replicates removed, stranded combined
test = t(ldply(unique(RNA_samples[which(RNA_samples$File %in% colnames(RNA_refseq_exon) & RNA_samples$Strand != "non"),]$Sample),function(x) apply(RNA_refseq_exon[,as.vector(RNA_samples[which(RNA_samples$Sample == x & (RNA_samples$Replicate == "r1a" | is.na(RNA_samples$Replicate))),]$File)],1,sum)))
colnames(test) = unique(RNA_samples[which(RNA_samples$File %in% colnames(RNA_refseq_exon) & RNA_samples$Strand != "non"),]$Sample)
RNA_refseq_exon_agnostic = cbind(RNA_refseq_exon[,1:6],test,RNA_refseq_exon[,105:108])
colnames(RNA_refseq_exon_agnostic)[56:59] = c("E028","E037","E038","E062")
rm(test)

# Number of samples exon is expressed >1 RPKM
RNA_refseq_exon_agnostic$Expressed_samples = apply(RNA_refseq_exon_agnostic[,8:59],1,function(x) sum(x > 1))

# Max expression per exon
RNA_refseq_exon_agnostic$Max_expression = apply(RNA_refseq_exon_agnostic[,8:59],1,max)

save(list=c("RNA_TE","RNA_TE_agnostic","RNA_refseq_exon","RNA_refseq_exon_agnostic"),file="R_datasets/rna.RData")