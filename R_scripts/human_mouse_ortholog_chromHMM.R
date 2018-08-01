# Human chromHMM state for hg19 TEs hypomethylated in mouse/human
hg19_orthologs_hypo_chromHMM = read.table("Mouse/hg19_ortholog_hypo_chromHMM.txt",sep='\t')[,c(1:8,10)]
colnames(hg19_orthologs_hypo_chromHMM) = c(hg19_coordinates,"Human_sample","Human_State")

# Mouse chromHMM state for mm10 TEs hypomethylated in mouse/human
mm10_orthologs_hypo_chromHMM = read.table("Mouse/mm10_ortholog_hypo_chromHMM.txt",sep='\t')[,c(1:8,10)]
colnames(mm10_orthologs_hypo_chromHMM) = c(mm10_coordinates,"Mouse_State","Mouse_sample_chromHMM")

# Combine chromHMM and WGBS
human_mouse_orthologs_WGBS_hypo_chromHMM = merge(human_mouse_orthologs_WGBS_hypo,hg19_orthologs_hypo_chromHMM,by=c(hg19_coordinates,"Human_sample"),all.x=TRUE)
human_mouse_orthologs_WGBS_hypo_chromHMM = merge(human_mouse_orthologs_WGBS_hypo_chromHMM,mm10_orthologs_hypo_chromHMM,by=c(mm10_coordinates,"Mouse_sample_chromHMM"),all.x=TRUE)
human_mouse_orthologs_WGBS_hypo_chromHMM$Mouse_State = factor(human_mouse_orthologs_WGBS_hypo_chromHMM$Mouse_State)
human_mouse_orthologs_WGBS_hypo_chromHMM$Human_State = factor(human_mouse_orthologs_WGBS_hypo_chromHMM$Human_State,chromHMM_states)
human_mouse_orthologs_WGBS_hypo_chromHMM = human_mouse_orthologs_WGBS_hypo_chromHMM[which(!is.na(human_mouse_orthologs_WGBS_hypo_chromHMM$Mouse_sample_chromHMM)),]