# Get samples where subfamilies are enriched in state
# 5/30/2017

# Samples where candidate subfamilies are enriched in state	 
#TE_landscape/enrichment/candidate_[state]_enrich.txt [5 files]	
write.table(subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "1_TssA" & subfamily_state_sample_filter$Subfamily %in% candidate_1TssA),c(1,5)],row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t',file="candidate_1TssA_enrich.txt")
write.table(subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "2_TssAFlnk" & subfamily_state_sample_filter$Subfamily %in% candidate_2TssAFlnk),c(1,5)],row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t',file="candidate_2TssAFlnk_enrich.txt")
write.table(subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "3_TxFlnk" & subfamily_state_sample_filter$Subfamily %in% candidate_3TxFlnk),c(1,5)],row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t',file="candidate_3TxFlnk_enrich.txt")
write.table(subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "6_EnhG" & subfamily_state_sample_filter$Subfamily %in% candidate_6EnhG),c(1,5)],row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t',file="candidate_6EnhG_enrich.txt")
write.table(subfamily_state_sample_filter[which(subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "7_Enh" & subfamily_state_sample_filter$Subfamily %in% candidate_7Enh),c(1,5)],row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t',file="candidate_7Enh_enrich.txt")