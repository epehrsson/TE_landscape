# Load RNA matrices (TEs and exons)
# See 3/1/2017, 3/2/2017, 3/3/2017, 3/7/2017, 6/8/2017, 6/11/2017, 6/14/2017, 8/1/2017

# Normalization factors per sample
RNA_metadata = read.table("raw_data/RNAseq/all.EGID.N.readlength",sep='\t',header=TRUE)

# RNA-seq sample metadata
RNA_samples = as.data.frame(read.table("sample_lists/RNA_samples_agnostic.txt"))
colnames(RNA_samples)[1] = "Sample"
RNA_samples = merge(RNA_samples,RNA_metadata,by.x="Sample",by.y="EGID")
RNA_samples$Factor = (1000/RNA_samples$read_length)/(RNA_samples$Total_number_of_reads_for_RPKM/1000000)
rm(RNA_metadata)

print("Loaded metadata")

# TEs
# Average expression per TE
RNA_TE = read.table("RNAseq/rmsk_TEother_average.txt",sep='\t')
colnames(RNA_TE) = c("chromosome","start","stop","subfamily","class","family","strand",as.vector(read.table("sample_lists/RNA_samples_agnostic.txt")$V1))
RNA_TE[,8:63] = t(t(RNA_TE[,8:63])*RNA_samples[match(colnames(RNA_TE[,8:63]),RNA_samples$Sample),]$Factor)

print("Loaded TEs")

# Number of samples TE is expressed >1 RPKM 
RNA_TE$Expressed_samples = apply(RNA_TE[,8:63],1,function(x) sum(x > 1))

# Number of samples TE is expressed >1 RPKM, no cancer cell lines/IMR90
RNA_TE$Expressed_samples_noCancer = apply(RNA_TE[,8:63],1,function(x) sum(as.numeric(x[which(metadata[match(colnames(RNA_TE)[8:63],metadata$Sample),]$Exclude == "Include")] > 1)))

# Max expression per TE
RNA_TE$Max_expression = apply(RNA_TE[,8:63],1,max)

# Updating class
RNA_TE$class_update = convert_class(RNA_TE$class)

print("Added columns")

# Write matrix for compare marks
write.table(melt(RNA_TE[,1:63],id.vars=TE_coordinates),file="RNAseq/rmsk_TE_rpkm.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')

print("Wrote matrix")

# Exons
# Average expression per Refseq exon
RNA_refseq_exon = read.table("RNAseq/refseq_exons_average.txt",sep='\t')
colnames(RNA_refseq_exon) = c("chromosome","start","stop","strand",as.vector(read.table("sample_lists/RNA_samples_agnostic.txt")$V1))
RNA_refseq_exon[,5:60] = t(t(RNA_refseq_exon[,5:60])*RNA_samples[match(colnames(RNA_refseq_exon[,5:60]),RNA_samples$Sample),]$Factor)

print("Loaded exons")

# Number of samples exon is expressed >1 RPKM
RNA_refseq_exon$Expressed_samples = apply(RNA_refseq_exon[,5:60],1,function(x) sum(x > 1))

# Max expression per exon
RNA_refseq_exon$Max_expression = apply(RNA_refseq_exon[,5:60],1,max)

print("Added columns")

save(list=c("RNA_TE","RNA_refseq_exon"),file="R_datasets/rna.RData")