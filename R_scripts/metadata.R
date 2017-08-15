# Metadata
# See 9/18/2016, 9/19/2016, 11/9/2016, 2/15/2017

# Metadata for human chromHMM samples

# Plus additional categories
EID_add = read.table("EID_metadata_addl.txt",header=TRUE,sep="\t")
colnames(EID_add)[1] = "Sample"
EID_metadata = merge(EID_metadata,EID_add,by=c("Sample"))

# Plus germlayer
EID_germ = read.table("EID_germ.txt",sep='\t')
colnames(EID_germ) = c("Sample","Germline")
EID_metadata = merge(EID_metadata,EID_germ,by=c("Sample"))

# Metadata for WGBS samples
EID_metadata_meth = EID_metadata[which(EID_metadata$Tissue %in% colnames(TE_meth_subfamily_num)),]

# Adding new categories
EID_metadata_meth = merge(EID_metadata_meth,EID_metadata[,c(1,9:10)],by=c("Sample"))
EID_metadata_meth = EID_metadata[match(EID_metadata_meth$Sample,EID_metadata$Sample),]
