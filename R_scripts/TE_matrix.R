# Create matrix of TEs with age, mappability, and feature overlaps
# See 4/19/2016, 4/25/2016, 8/24/2016, 8/25/2016, 9/20/2016, 9/21/2016, 9/27/2016, 9/28/2016, 11/4/2016, 11/5/2016, 11/7/2016, 11/18/2016, 12/16/2016, 1/31/2017, 2/1/2017, 2/3/2017, 2/6/2017, 2/9/2017, 2/10/2017, 2/25/2017, 2/27/2017, 2/28/2017, 3/5/2017, 3/8/2017, 5/14/2017, 5/15/2017, 5/16/2017, 5/17/2017, 6/7/2017, 6/14/2017, 6/15/2017, 7/21/2017, 7/24/2017, 8/1/2017, 8/2/2017

library(plyr)
library(reshape2)

# RepeatMasker file restricted to TEs
rmsk_TE = read.table(file="features/TEs/rmsk_TEother.txt",sep='\t')
colnames(rmsk_TE) = c("chromosome","start","stop","subfamily","class","family","strand")
rmsk_TE$Length = rmsk_TE$stop-rmsk_TE$start

rmsk_TE$class_update = rmsk_TE$class
rmsk_TE$class_update = factor(rmsk_TE$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
rmsk_TE[which(rmsk_TE$class == "Other"),]$class_update = "SVA"
rmsk_TE[which(rmsk_TE$class %in% c("DNA?","LINE?","SINE?","LTR?","Unknown","Unknown?","RC")),]$class_update = "Other"

# Add mappability per TE
rmsk_TE_map36 = rbind(read.table("mappability/rmsk_TE_mappability_36mer.txt",sep='\t'),read.table("mappability/rmsk_other_mappability_36mer.txt",sep='\t'))
colnames(rmsk_TE_map36) = c("chromosome","start","stop","subfamily","mappability")
rmsk_TE = merge(rmsk_TE,rmsk_TE_map36,by=c("chromosome","start","stop","subfamily"))
rm(rmsk_TE_map36)

# Add age (JC evolutionary distance) per TE
rmsk_TE_age = rbind(read.table("age/rmsk_TE_JCage.txt",sep='\t'),rmsk_other_age = read.table("age/rmsk_other_JCage.txt",sep='\t'))
colnames(rmsk_TE_age) = c("chromosome","start","stop","subfamily","family","class","strand","substitutions","JC_distance")
rmsk_TE = merge(rmsk_TE,rmsk_TE_age,by=c("chromosome","start","stop","subfamily","family","class","strand"))
rm(rmsk_TE_age)

# Overlap between TEs and Refseq features
feature_files = list.files(path="features/intersect_features/",pattern="rmsk_TEother_",full.names=TRUE)
features = lapply(feature_files,function(x) read.table(x,sep='\t'))
features = lapply(features, setNames, nm =c("chromosome","start","stop","subfamily","class","family","strand"))
feature_files = gsub("features/intersect_features//rmsk_TEother_","",feature_files)
feature_files = gsub("refseq_","",feature_files)
feature_files = gsub("_merge.txt","",feature_files)
feature_files = gsub(".txt","",feature_files)
for (i in 1:length(feature_files)){
  colnames(features[[i]])[8] = feature_files[i]
}
features_merge = reshape::merge_recurse(features,by=c("chromosome","start","stop","subfamily","family","class","strand"),all.x=TRUE)

# Convert to matrix
rmsk_TE = merge(rmsk_TE,features_merge,by=c("chromosome","start","stop","subfamily","family","class","strand"),all.x=TRUE)
rm(list=c("feature_files","features","features_merge"))

# Add TEs overlapping the hg19 blacklist
TE_blacklist = read.table("features/intersect_features/TE_blacklist.bed",sep='\t')
colnames(TE_blacklist) = c("chromosome","start","stop","subfamily","class","family","strand","chromosome_blacklist","start_blacklist","stop_blacklist","blacklist")
TE_blacklist = aggregate(data=TE_blacklist,blacklist~chromosome+start+stop+subfamily+class+family+strand,sum)
rmsk_TE = merge(rmsk_TE,TE_blacklist,by=c("chromosome","start","stop","subfamily","class","family","strand"),all.x=TRUE)
rm(TE_blacklist)

save(rmsk_TE,file="R_scripts/rmsk_TE.RData")
