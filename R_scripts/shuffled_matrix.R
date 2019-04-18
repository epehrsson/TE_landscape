# Load shuffled TE matrices

print("Loading shuffled matrices")
rmsk_TE_shuffled = lapply(seq(1,10,1),function(x) read.table(paste("features/shuffled_TEs/rmsk_TE_shuffle_",x,".txt",sep=""),sep='\t',
                                                             col.names=TE_coordinates[c(1:4,6,5,7)]))

# Add shuffled TE features
print("Adding features")

for (i in 1:10){
  print(i)
  print("Loading features")
  feature_files = list.files(path="features/shuffled_TEs/intersect_features/",pattern=paste("rmsk_TE_shuffle_",i,"_",sep=""),full.names=TRUE)
  features = lapply(feature_files,function(x) read.table(x,sep='\t',col.names=c(TE_coordinates[c(1:4,6,5,7)],"X")))
  feature_files = gsub(paste("features/shuffled_TEs/intersect_features//rmsk_TE_shuffle_",i,"_",sep=""),"",feature_files)
  feature_files = gsub("refseq_","",feature_files)
  feature_files = gsub("_merge.txt","",feature_files)
  feature_files = gsub(".txt","",feature_files)
  for (j in 1:length(feature_files)){
    colnames(features[[j]])[8] = feature_files[j]
  }
  
  print("Merging features")
  features_merge = reshape::merge_recurse(features,by=TE_coordinates,all.x=TRUE)
  
  print("Adding feature overlap")
  rmsk_TE_shuffled[[i]] = merge(rmsk_TE_shuffled[[i]],features_merge,by=TE_coordinates,all.x=TRUE)
  
  rm(list=c("feature_files","features","features_merge"))
}

save(rmsk_TE_shuffled,file="R_datasets/shuffled_TEs.RData")