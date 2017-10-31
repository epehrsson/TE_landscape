# Subfamily enrichment candidate subfamiles
# See 7/11/2016, 8/29/2016, 9/29/2016, 11/22/2016, 11/27/2016, 5/23/2017, 5/24/2017, 5/29/2017, 5/30/2017, 6/15/2017

# Statistics for candidate subfamilies
test = merge(merge(aggregate(data=subfamily_state_sample_filter[which(subfamily_state_sample_filter$Subfamily %in% candidate_1TssA & subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "1_TssA"),],Members~Subfamily,min),aggregate(data=subfamily_state_sample_filter[which(subfamily_state_sample_filter$Subfamily %in% candidate_1TssA & subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "1_TssA"),],Members~Subfamily,max),by=c("Subfamily")),aggregate(data=subfamily_state_sample_filter[which(subfamily_state_sample_filter$Subfamily %in% candidate_1TssA & subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "1_TssA"),],Members~Subfamily,median),by=c("Subfamily"))
colnames(test) = c("Subfamily","Min","Max","Median")
test = merge(test,aggregate(data=potential_TEother_state[which(potential_TEother_state$subfamily %in% candidate_1TssA),],X1_TssA~subfamily,function(x) sum(x > 0)),by.x=c("Subfamily"),by.y=c("subfamily"))
test = merge(test,subfamily_state_sample_counts[which(subfamily_state_sample_counts$State == "1_TssA"),c(3,5)],by=c("Subfamily"))
colnames(test)[5:6] = c("Members_ever","Sample_enriched")
test = merge(test,rmsk_TEother_stats_subfamily[,3:4],by=c("Subfamily"))
colnames(test)[7] = "Members"
test = merge(test,aggregate(data=potential_TEother_state[which(potential_TEother_state$subfamily %in% candidate_1TssA),],X1_TssA~subfamily,function(x) sum(x > 5)),by.x=c("Subfamily"),by.y=c("subfamily"))

# Statistics for candidate subfamilies
test = merge(merge(aggregate(data=subfamily_state_sample_filter[which(subfamily_state_sample_filter$Subfamily %in% candidate_7Enh & subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "7_Enh"),],Members~Subfamily,min),aggregate(data=subfamily_state_sample_filter[which(subfamily_state_sample_filter$Subfamily %in% candidate_7Enh & subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "7_Enh"),],Members~Subfamily,max),by=c("Subfamily")),aggregate(data=subfamily_state_sample_filter[which(subfamily_state_sample_filter$Subfamily %in% candidate_7Enh & subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "7_Enh"),],Members~Subfamily,median),by=c("Subfamily"))
colnames(test) = c("Subfamily","Min","Max","Median")
test = merge(test,aggregate(data=potential_TEother_state[which(potential_TEother_state$subfamily %in% candidate_7Enh),],X7_Enh~subfamily,function(x) sum(x > 0)),by.x=c("Subfamily"),by.y=c("subfamily"))
test = merge(test,subfamily_state_sample_counts[which(subfamily_state_sample_counts$State == "7_Enh"),c(3,5)],by=c("Subfamily"))
colnames(test)[5:6] = c("Members_ever","Sample_enriched")
test = merge(test,rmsk_TEother_stats_subfamily[,3:4],by=c("Subfamily"))
colnames(test)[7] = "Members"
test = merge(test,aggregate(data=potential_TEother_state[which(potential_TEother_state$subfamily %in% candidate_7Enh),],X7_Enh~subfamily,function(x) sum(x > 5)),by.x=c("Subfamily"),by.y=c("subfamily"))

# Write enriched subfamily coordinates and enriched samples
# TEs ever in state, by subfamily	 
#TE_landscape/enrichment/*_1TssA.bed [40 files]	
lapply(candidate_1TssA,function(x) write.table(potential_TEother_state[which(potential_TEother_state$subfamily == x & potential_TEother_state$X1_TssA > 0),c(1:4,7)],sep='\t',row.names=FALSE,quote=FALSE,file=paste("enrichment/",x,"_1TssA.bed",sep="")))
#TE_landscape/enrichment/*_2TssAFlnk.bed [93 files]		 
lapply(candidate_2TssAFlnk,function(x) write.table(potential_TEother_state[which(potential_TEother_state$subfamily == x & potential_TEother_state$X2_TssAFlnk > 0),c(1:4,9,7)],sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,file=paste("enrichment/",x,"_2TssAFlnk.bed",sep="")))
#TE_landscape/enrichment/*_3TxFlnk.bed [14 files]		 
lapply(candidate_3TxFlnk,function(x) write.table(potential_TEother_state[which(potential_TEother_state$subfamily == x & potential_TEother_state$X3_TxFlnk > 0),c(1:4,10,7)],sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,file=paste("enrichment/",x,"_3TxFlnk.bed",sep="")))
#TE_landscape/enrichment/*_6EnhG.bed [38 files]		 
lapply(candidate_6EnhG,function(x) write.table(potential_TEother_state[which(potential_TEother_state$subfamily == x & potential_TEother_state$X6_EnhG > 0),c(1:4,13,7)],sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,file=paste("enrichment/",x,"_6EnhG.bed",sep="")))
#TE_landscape/enrichment/*_7Enh.bed [166 files]		 
lapply(candidate_7Enh,function(x) write.table(potential_TEother_state[which(potential_TEother_state$subfamily == x & potential_TEother_state$X7_Enh > 0),c(1:3,4,14,7)],sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,file=paste("enrichment/",x,"_7Enh.bed",sep="")))

# TEs never in state, by subfamily	 
#TE_landscape/enrichment/*_no1TssA.bed [40 files]	
lapply(candidate_1TssA,function(x) write.table(potential_TEother_state[which(potential_TEother_state$subfamily == x & potential_TEother_state$X1_TssA == 0),c(1:4,7)],sep='\t',row.names=FALSE,quote=FALSE,file=paste("enrichment/",x,"_no1TssA.bed",sep="")))
#TE_landscape/enrichment/*_no2TssAFlnk.bed [93 files]		 
lapply(candidate_2TssAFlnk,function(x) write.table(potential_TEother_state[which(potential_TEother_state$subfamily == x & potential_TEother_state$X2_TssAFlnk == 0),c(1:4,9,7)],sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,file=paste("enrichment/",x,"_no2TssAFlnk.bed",sep="")))
#TE_landscape/enrichment/*_no3TxFlnk.bed [14 files]		
lapply(candidate_3TxFlnk,function(x) write.table(potential_TEother_state[which(potential_TEother_state$subfamily == x & potential_TEother_state$X3_TxFlnk == 0),c(1:4,10,7)],sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,file=paste("enrichment/",x,"_no3TxFlnk.bed",sep="")))
#TE_landscape/enrichment/*_no6EnhG.bed [38 files]		 
lapply(candidate_6EnhG,function(x) write.table(potential_TEother_state[which(potential_TEother_state$subfamily == x & potential_TEother_state$X6_EnhG == 0),c(1:4,13,7)],sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,file=paste("enrichment/",x,"_no6EnhG.bed",sep="")))
#TE_landscape/enrichment/*_no7Enh.bed [166 files]		 
lapply(candidate_7Enh,function(x) write.table(potential_TEother_state[which(potential_TEother_state$subfamily == x & potential_TEother_state$X7_Enh == 0),c(1:3,4,14,7)],sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,file=paste("enrichment/",x,"_no7Enh.bed",sep="")))

# Samples where candidate subfamilies are enriched in state	 
#TE_landscape/enrichment/candidate_[state]_enrich.txt [5 files]	
write.table(subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "1_TssA" & subfamily_state_sample_filter$Subfamily %in% candidate_1TssA),c(1,5)],row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t',file="candidate_1TssA_enrich.txt")
write.table(subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "2_TssAFlnk" & subfamily_state_sample_filter$Subfamily %in% candidate_2TssAFlnk),c(1,5)],row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t',file="candidate_2TssAFlnk_enrich.txt")
write.table(subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "3_TxFlnk" & subfamily_state_sample_filter$Subfamily %in% candidate_3TxFlnk),c(1,5)],row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t',file="candidate_3TxFlnk_enrich.txt")
write.table(subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "6_EnhG" & subfamily_state_sample_filter$Subfamily %in% candidate_6EnhG),c(1,5)],row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t',file="candidate_6EnhG_enrich.txt")
write.table(subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "7_Enh" & subfamily_state_sample_filter$Subfamily %in% candidate_7Enh),c(1,5)],row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t',file="candidate_7Enh_enrich.txt")