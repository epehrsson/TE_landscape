# chromHMM
## Human chromHMM state for hg19 orthologous TEs
hg19_orthologs_chromHMM = read.table("Mouse/chromHMM/hg19_chromHMM_TE.txt",sep='\t')[,c(1:8,10)]
colnames(hg19_orthologs_chromHMM) = c(hg19_coordinates,"Human_sample","Human_state_chromHMM")

## Mouse chromHMM state for mm10 orthologous TEs
mm10_orthologs_chromHMM = read.table("Mouse/chromHMM/mm10_chromHMM_TE_sorted.txt",sep='\t')[,1:9]
colnames(mm10_orthologs_chromHMM) = c(mm10_coordinates,"Mouse_state_chromHMM","Mouse_sample_chromHMM")

## Orthologous TEs with chromHMM state, paired samples only
human_mouse_orthologs_chromHMM_paired = apply(human_mouse_samples,1,function(x) merge(merge(human_mouse_orthologs_mm10[,c(hg19_coordinates,mm10_coordinates,"class_update","JC_distance")],
                                                                                            hg19_orthologs_chromHMM[which(hg19_orthologs_chromHMM$Human_sample == x[1]),],by=hg19_coordinates),
                                                                                      mm10_orthologs_chromHMM[which(mm10_orthologs_chromHMM$Mouse_sample_chromHMM == x[3]),],by=mm10_coordinates))
human_mouse_orthologs_chromHMM_paired = ldply(human_mouse_orthologs_chromHMM_paired)
human_mouse_orthologs_chromHMM_paired$Human_state_chromHMM = factor(human_mouse_orthologs_chromHMM_paired$Human_state_chromHMM,chromHMM_states)
human_mouse_orthologs_chromHMM_paired$Mouse_state_chromHMM = factor(human_mouse_orthologs_chromHMM_paired$Mouse_state_chromHMM,mouse_chromHMM_states)
