# chromHMM
## Human chromHMM state for hg19 orthologous TEs
hg19_orthologs_chromHMM = read.table("Mouse/chromHMM/hg19_chromHMM_TE.txt",sep='\t')[,c(1:8,10)]
colnames(hg19_orthologs_chromHMM) = c(hg19_coordinates,"Human_sample","Human_state_chromHMM")
hg19_orthologs_chromHMM = merge(hg19_orthologs_chromHMM,human_mouse_samples[,c("Human_sample",sample_categories)],by="Human_sample")

## Mouse chromHMM state for mm10 orthologous TEs
mm10_orthologs_chromHMM = read.table("Mouse/chromHMM/mm10_chromHMM_TE_sorted.txt",sep='\t')[,1:9]
colnames(mm10_orthologs_chromHMM) = c(mm10_coordinates,"Mouse_state_chromHMM","Mouse_sample_chromHMM")
mm10_orthologs_chromHMM = merge(mm10_orthologs_chromHMM,human_mouse_samples[,c("Mouse_sample_chromHMM",sample_categories)],by="Mouse_sample_chromHMM")

## Orthologous TEs with chromHMM state, paired samples only
human_mouse_orthologs_chromHMM_paired = apply(human_mouse_samples,1,function(x) merge(merge(human_mouse_orthologs_mm10[,c(hg19_coordinates,mm10_coordinates,"class_update","JC_distance")],
                                                                                            hg19_orthologs_chromHMM[which(hg19_orthologs_chromHMM$Human_sample == x[1]),],by=hg19_coordinates),
                                                                                      mm10_orthologs_chromHMM[which(mm10_orthologs_chromHMM$Mouse_sample_chromHMM == x[3]),],by=mm10_coordinates))
human_mouse_orthologs_chromHMM_paired = ldply(human_mouse_orthologs_chromHMM_paired)
human_mouse_orthologs_chromHMM_paired$Human_state_chromHMM = factor(human_mouse_orthologs_chromHMM_paired$Human_state_chromHMM,chromHMM_states)
human_mouse_orthologs_chromHMM_paired$Mouse_state_chromHMM = factor(human_mouse_orthologs_chromHMM_paired$Mouse_state_chromHMM,mouse_chromHMM_states)

## Tissue-specific TEs, active regulatory states
test1 = ddply(hg19_orthologs_chromHMM[which(hg19_orthologs_chromHMM$Human_state_chromHMM %in% chromHMM_states[c(1:3,6:7)]),],
              .(human_chr_hg19,human_start_hg19,human_stop_hg19,human_subfamily,human_class,human_family,human_strand_hg19,Anatomy),
              summarise,Human_samples=length(unique(Human_sample)))
test1 = dcast(test1,human_chr_hg19+human_start_hg19+human_stop_hg19+human_subfamily+human_class+human_family+human_strand_hg19~Anatomy,value.var="Human_samples")

test2 = ddply(mm10_orthologs_chromHMM[which(mm10_orthologs_chromHMM$Mouse_state_chromHMM %in% mouse_chromHMM_states[c(1:3,6:8)]),],
              .(mouse_chr_mm10,mouse_start_mm10,mouse_stop_mm10,mouse_subfamily,mouse_class,mouse_family,mouse_strand_mm10,Anatomy),
              summarise,Mouse_samples=length(unique(Mouse_sample_chromHMM)))
test2 = dcast(test2,mouse_chr_mm10+mouse_start_mm10+mouse_stop_mm10+mouse_subfamily+mouse_class+mouse_family+mouse_strand_mm10~Anatomy,value.var="Mouse_samples")

hg19_mm10_chromHMM_active = merge(test1,human_mouse_orthologs_mm10[,c(hg19_coordinates,mm10_coordinates,"JC_distance")],by=hg19_coordinates,all.y=TRUE)
hg19_mm10_chromHMM_active = merge(hg19_mm10_chromHMM_active,test2,by=mm10_coordinates,all.x=TRUE)
rm(list=c("test1","test2"))

hg19_mm10_chromHMM_active[is.na(hg19_mm10_chromHMM_active)] = 0
hg19_mm10_chromHMM_active$GI_INTESTINE.x = hg19_mm10_chromHMM_active$GI_INTESTINE.x + hg19_mm10_chromHMM_active$GI_COLON.x
hg19_mm10_chromHMM_active$GI_INTESTINE.y = hg19_mm10_chromHMM_active$GI_INTESTINE.y + hg19_mm10_chromHMM_active$GI_COLON.y
hg19_mm10_chromHMM_active = hg19_mm10_chromHMM_active[,c(1:14,23,15,17:22,24,26:31)]
hg19_mm10_chromHMM_active$Human_samples = apply(hg19_mm10_chromHMM_active,1,function(x) sum(as.numeric(x[16:22])))
hg19_mm10_chromHMM_active$Mouse_samples = apply(hg19_mm10_chromHMM_active,1,function(x) sum(as.numeric(x[23:29])))

## Tissue-specific TEs, promoters
test1 = ddply(hg19_orthologs_chromHMM[which(hg19_orthologs_chromHMM$Human_state_chromHMM == "1_TssA"),],
              .(human_chr_hg19,human_start_hg19,human_stop_hg19,human_subfamily,human_class,human_family,human_strand_hg19,Anatomy),
              summarise,Human_samples=length(unique(Human_sample)))
test1 = dcast(test1,human_chr_hg19+human_start_hg19+human_stop_hg19+human_subfamily+human_class+human_family+human_strand_hg19~Anatomy,value.var="Human_samples")

test2 = ddply(mm10_orthologs_chromHMM[which(mm10_orthologs_chromHMM$Mouse_state_chromHMM == "TssA"),],
              .(mouse_chr_mm10,mouse_start_mm10,mouse_stop_mm10,mouse_subfamily,mouse_class,mouse_family,mouse_strand_mm10,Anatomy),
              summarise,Mouse_samples=length(unique(Mouse_sample_chromHMM)))
test2 = dcast(test2,mouse_chr_mm10+mouse_start_mm10+mouse_stop_mm10+mouse_subfamily+mouse_class+mouse_family+mouse_strand_mm10~Anatomy,value.var="Mouse_samples")

hg19_mm10_chromHMM_promoter = merge(test1,human_mouse_orthologs_mm10[,c(hg19_coordinates,mm10_coordinates,"JC_distance")],by=hg19_coordinates,all.y=TRUE)
hg19_mm10_chromHMM_promoter = merge(hg19_mm10_chromHMM_promoter,test2,by=mm10_coordinates,all.x=TRUE)
rm(list=c("test1","test2"))

hg19_mm10_chromHMM_promoter[is.na(hg19_mm10_chromHMM_promoter)] = 0
hg19_mm10_chromHMM_promoter$GI_INTESTINE.x = hg19_mm10_chromHMM_promoter$GI_INTESTINE.x + hg19_mm10_chromHMM_promoter$GI_COLON.x
hg19_mm10_chromHMM_promoter$GI_INTESTINE.y = hg19_mm10_chromHMM_promoter$GI_INTESTINE.y + hg19_mm10_chromHMM_promoter$GI_COLON.y
hg19_mm10_chromHMM_promoter = hg19_mm10_chromHMM_promoter[,c(1:14,23,15,17:22,24,26:31)]
hg19_mm10_chromHMM_promoter$Human_samples = apply(hg19_mm10_chromHMM_promoter,1,function(x) sum(as.numeric(x[16:22])))
hg19_mm10_chromHMM_promoter$Mouse_samples = apply(hg19_mm10_chromHMM_promoter,1,function(x) sum(as.numeric(x[23:29])))
