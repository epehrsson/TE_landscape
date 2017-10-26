# Mouse-human candidates
# 5/30/2017, 7/26/2017, 7/31/2017, 10/3/2017

# hg19 TEs hypomethlated in human/mouse also in promoter state	
#TE_landscape/Mouse/human_mouse_hypo_1TssA.bed	
write.table(unique(human_mouse_orthologs_mm10_hypo[which(human_mouse_orthologs_mm10_hypo$Human_state == "1_TssA" & human_mouse_orthologs_mm10_hypo$Mouse_State == "1"),1:7]),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t',file="Mouse/human_mouse_hypo_1TssA.bed")

# Overlap with GENCODE promoters	
#TE_landscape/Mouse/human_mouse_hypo_1TssA_Gencode.txt	
bedtools intersect -wo -a human_mouse_hypo_1TssA.bed -b ~/genic_features/Gencodev19_promoter.txt > human_mouse_hypo_1TssA_Gencode.txt

# Update 10/3/2017
## Overlap with GENCODE promoters 
bedtools intersect -wo -a hg19_orthologs_hypo_1TssA.bed -b ~/genic_features/Gencode/Gencodev19_promoter.txt > hg19_orthologs_hypo_1TssA_Gencode.txt
