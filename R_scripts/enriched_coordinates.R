# Write enriched subfamily coordinates
# 5/23/2017, 5/29/2017, 5/30/2017

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
