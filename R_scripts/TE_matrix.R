# Creates a matrix of individual hg19 TEs ("rmsk_TE")
# Where rows are TEs (n=4,430,788) and columns are TE characteristics such as age, mappability, and feature overlap

# Load RepeatMasker file, restricted to TEs
# Columns are TE chromosome, start, stop and strand, plus subfamily, family, and class
rmsk_TE = read.table(file="features/TEs/rmsk_TEother.txt",sep='\t',col.names=TE_coordinates[c(1:4,6,5,7)])

# Add length of TE
rmsk_TE$Length = rmsk_TE$stop-rmsk_TE$start

# Update class assignment
rmsk_TE$class_update = convert_class(rmsk_TE$class)

# Add mappability per TE
rmsk_TE_map36 = rbind(read.table("mappability/rmsk_TE_mappability_36mer.txt",sep='\t'),read.table("mappability/rmsk_other_mappability_36mer.txt",sep='\t'))
colnames(rmsk_TE_map36) = c("chromosome","start","stop","subfamily","mappability")
rmsk_TE = merge(rmsk_TE,rmsk_TE_map36,by=c("chromosome","start","stop","subfamily"))
rm(rmsk_TE_map36)

# Add age (Jukes-Cantor evolutionary distance) per TE
rmsk_TE_age = rbind(read.table("age/rmsk_TE_JCage.txt",sep='\t'),rmsk_other_age = read.table("age/rmsk_other_JCage.txt",sep='\t'))
colnames(rmsk_TE_age) = c("chromosome","start","stop","subfamily","family","class","strand","substitutions","JC_distance")
rmsk_TE = merge(rmsk_TE,rmsk_TE_age,by=TE_coordinates)
rm(rmsk_TE_age)

# Add number of CpGs per TE
TE_CpG_count = read.table("WGBS/TE_CpG_count.txt",sep='\t',col.names=c(TE_coordinates[c(1:4,6,5,7)],"CpGs"))
TE_CpG_count$CpGs = TE_CpG_count$CpGs/2
rmsk_TE = merge(rmsk_TE,TE_CpG_count,all.x=TRUE,by=TE_coordinates)
rmsk_TE[is.na(rmsk_TE)] = 0
rm(TE_CpG_count)

# Add CpG density
rmsk_TE$CpGs_per_length = rmsk_TE$CpGs/rmsk_TE$Length

# Add length of overlap between TEs and RefSeq gene features
feature_files = list.files(path="features/intersect_features/",pattern="rmsk_TEother_",full.names=TRUE)
features = lapply(feature_files,function(x) read.table(x,sep='\t'))
features = lapply(features, setNames, nm =TE_coordinates[c(1:4,6,5,7)])
feature_files = gsub("features/intersect_features//rmsk_TEother_","",feature_files)
feature_files = gsub("refseq_","",feature_files)
feature_files = gsub("_merge.txt","",feature_files)
feature_files = gsub(".txt","",feature_files)
for (i in 1:length(feature_files)){
  colnames(features[[i]])[8] = feature_files[i]
}
## Combine matrices
features_merge = reshape::merge_recurse(features,by=TE_coordinates,all.x=TRUE)

## Add to TE matrix
rmsk_TE = merge(rmsk_TE,features_merge,by=TE_coordinates,all.x=TRUE)
rm(list=c("feature_files","features","features_merge"))

# Add length of overlap with hg19 blacklist
TE_blacklist = read.table("features/intersect_features/TE_blacklist.bed",sep='\t',
                          col.names=c(TE_coordinates[c(1:4,6,5,7)],"chromosome_blacklist","start_blacklist","stop_blacklist","blacklist"))
TE_blacklist = aggregate(data=TE_blacklist,blacklist~chromosome+start+stop+subfamily+class+family+strand,sum)
rmsk_TE = merge(rmsk_TE,TE_blacklist,by=TE_coordinates,all.x=TRUE)
rm(TE_blacklist)

# Add overlap with VISTA enhancers
vista_TE = read.table("features/intersect_features/vista_enhancers_TE.txt",sep='\t',
                      col.names=c(TE_coordinates[c(1:4,6,5,7)],"Vista_enhancers"))
rmsk_TE = merge(rmsk_TE,vista_TE,by=TE_coordinates,all.x=TRUE)
rm(vista_TE)

# Save dataframe
save(rmsk_TE,file="R_datasets/rmsk_TE.RData")