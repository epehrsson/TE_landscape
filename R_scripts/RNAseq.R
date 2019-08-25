# Creates matrices with the expression level of each TE or exon (row) per sample (column)
# For all samples with strand-agnostic RNA-seq data (n=56)
## RNA_TE: TEs
## RNA_refseq_exon: RefSeq exons

# Load the total reads and read length for each sample, provided by Roadmap
RNA_metadata = read.table("raw_data/RNAseq/all.EGID.N.readlength",sep='\t',header=TRUE)

# Calculate the normalization factor by which read coverage is scaled to produce RPKM
RNA_samples = as.data.frame(read.table("sample_lists/RNA_samples_agnostic.txt"))
colnames(RNA_samples)[1] = "Sample"
RNA_samples = merge(RNA_samples,RNA_metadata,by.x="Sample",by.y="EGID")
RNA_samples$Factor = (1000/RNA_samples$read_length)/(RNA_samples$Total_number_of_reads_for_RPKM/1000000)
rm(RNA_metadata)

print("Loaded metadata")

# Average expression per TE
## Load RNA-seq read coverage per TE
RNA_TE = read.table("RNAseq/rmsk_TEother_average.txt",sep='\t')
colnames(RNA_TE) = c("chromosome","start","stop","subfamily","class","family","strand",as.vector(read.table("sample_lists/RNA_samples_agnostic.txt")$V1))

## Convert to RPKM
RNA_TE[,8:63] = t(t(RNA_TE[,8:63])*RNA_samples[match(colnames(RNA_TE[,8:63]),RNA_samples$Sample),]$Factor)

print("Loaded TEs")

## Number of samples the TE is expressed >1 RPKM 
RNA_TE$Expressed_samples = apply(RNA_TE[,8:63],1,function(x) sum(x > 1))

## Number of samples the TE is expressed >1 RPKM, excluding cancer cell lines/IMR90
RNA_TE$Expressed_samples_noCancer = apply(RNA_TE[,8:63],1,function(x) sum(as.numeric(x[which(metadata[match(colnames(RNA_TE)[8:63],metadata$Sample),]$Exclude == "Include")] > 1)))

## Max expression per TE
RNA_TE$Max_expression = apply(RNA_TE[,8:63],1,max)

## Update class assignments
RNA_TE$class_update = convert_class(RNA_TE$class)

print("Added columns")

# Write out matrix
write.table(melt(RNA_TE[,1:63],id.vars=TE_coordinates),file="RNAseq/rmsk_TE_rpkm.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')

print("Wrote matrix")

# Average expression per RefSeq exon
## Load RNA-seq read coverage per exon
RNA_refseq_exon = read.table("RNAseq/refseq_exons_average.txt",sep='\t')
colnames(RNA_refseq_exon) = c("chromosome","start","stop","strand",as.vector(read.table("sample_lists/RNA_samples_agnostic.txt")$V1))

## Convert to RPKM
RNA_refseq_exon[,5:60] = t(t(RNA_refseq_exon[,5:60])*RNA_samples[match(colnames(RNA_refseq_exon[,5:60]),RNA_samples$Sample),]$Factor)

print("Loaded exons")

## Number of samples the exon is expressed >1 RPKM
RNA_refseq_exon$Expressed_samples = apply(RNA_refseq_exon[,5:60],1,function(x) sum(x > 1))

## Max expression per exon
RNA_refseq_exon$Max_expression = apply(RNA_refseq_exon[,5:60],1,max)

print("Added columns")

# Save dataframes
save(list=c("RNA_TE","RNA_refseq_exon"),file="R_datasets/rna.RData")