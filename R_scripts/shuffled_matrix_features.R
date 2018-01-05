# Add shuffled TE features

load("R_datasets/shuffled_TEs_map.RData")
rmsk_TE_shuffled = rmsk_TE_shuffled[c(1,3:10,2)]

# Overlap between TEs and Refseq features
for (i in 1:10){
  print(i)
  feature_files = list.files(path="features/shuffled_TEs/intersect_features/",pattern=paste("rmsk_TE_shuffle_",i,"_",sep=""),full.names=TRUE)
  features = lapply(feature_files,function(x) read.table(x,sep='\t'))
  features = lapply(features, setNames, nm =c("chromosome","start","stop","subfamily","class","family","strand"))
  feature_files = gsub(paste("features/shuffled_TEs/intersect_features//rmsk_TE_shuffle_",i,"_",sep=""),"",feature_files)
  feature_files = gsub("refseq_","",feature_files)
  feature_files = gsub("_merge.txt","",feature_files)
  feature_files = gsub(".txt","",feature_files)
  for (j in 1:length(feature_files)){
    colnames(features[[j]])[8] = feature_files[j]
  }
  print("Loaded features")
  features_merge = reshape::merge_recurse(features,by=c("chromosome","start","stop","subfamily","family","class","strand"),all.x=TRUE)
  print("Merged features")
  rmsk_TE_shuffled[[i]] = merge(rmsk_TE_shuffled[[i]],features_merge,by=c("chromosome","start","stop","subfamily","family","class","strand"),all.x=TRUE)
  print("Added feature overlap")
  rm(list=c("feature_files","features","features_merge"))
}

save(rmsk_TE_shuffled,file="R_datasets/shuffled_TEs_feature.RData")