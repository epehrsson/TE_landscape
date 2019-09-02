# Creates a matrix of metadata ("metadata") for all consolidated Roadmap epigenomes
# Where row is epigenome/sample (n=127) and columns are metadata

# Categories from Roadmap consortium (for each EID, Group, Anatomy, and Type plus sample name and group color)
metadata = read.table("metadata/EID_metadata.txt",header=TRUE,sep='\t',comment.char = "")[,c(1,3:5,7:8)]
colnames(metadata)[1] = "Sample"

# Add Age and Cancer status
EID_add = read.table("metadata/EID_metadata_addl.txt",header=TRUE,sep="\t")
colnames(EID_add)[1] = "Sample"
metadata = merge(metadata,EID_add,by=c("Sample"))
rm(EID_add)

# Add Germlayer
EID_germ = read.table("metadata/EID_germ.txt",sep='\t',col.names=c("Sample","Germline"))
metadata = merge(metadata,EID_germ,by=c("Sample"))
rm(EID_germ)

# Add columns indicating whether sample has each data type
## WGBS
WGBS_samples = read.table("sample_lists/WGBS_samples.txt")
WGBS_samples$WGBS = rep("WGBS",dim(WGBS_samples)[2])
metadata = merge(metadata,WGBS_samples,by.x="Sample",by.y="V1",all.x=TRUE)
rm(WGBS_samples)

## DHS
DNase_samples = read.table("sample_lists/DNase_samples.txt")
DNase_samples$DNase = rep("DNase",dim(DNase_samples)[2])
metadata = merge(metadata,DNase_samples,by.x="Sample",by.y="V1",all.x=TRUE)
rm(DNase_samples)

## H3K27ac
H3K27ac_samples = read.table("sample_lists/H3K27ac_samples.txt")
H3K27ac_samples$H3K27ac = rep("H3K27ac",dim(H3K27ac_samples)[2])
metadata = merge(metadata,H3K27ac_samples,by.x="Sample",by.y="V1",all.x=TRUE)
rm(H3K27ac_samples)

## RNA-seq
RNA_samples = read.table("sample_lists/RNA_samples_agnostic.txt")
RNA_samples$RNA = rep("RNA",dim(RNA_samples)[2])
metadata = merge(metadata,RNA_samples,by.x="Sample",by.y="V1",all.x=TRUE)
rm(RNA_samples)

# Add column indicating whether sample includes chrY
metadata$chrY = rep("Yes",127)
metadata[which(metadata$Sample %in% c("E116","E117","E123","E124","E126","E127")),]$chrY = "No"

# Add column indicating whether sample is IMR90 or cancer cell line
metadata$Exclude = rep("Include",127)
metadata[which(metadata$Sample %in% c("E017","E114","E115","E117","E118","E123")),]$Exclude = "Exclude"