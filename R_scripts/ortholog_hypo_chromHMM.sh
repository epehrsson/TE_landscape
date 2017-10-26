# Human-mouse hypomethylated orthologs, chromHMM
# 5/29/2017, 5/30/2017, 7/26/2017, 7/31/2017, 10/3/2017

# hg19 TEs hypomethylated in human and mouse	
#TE_landscape/human_mouse_ortholog_hypo.txt	
write.table(as.matrix(unique(human_mouse_orthologs_mm10_hypo[,c(2:8,16)])),file="human_mouse_ortholog_hypo.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

# Human chromHMM state for hg19 TEs hypomethylated in mouse/human	
#TE_landscape/Mouse/human_mouse_ortholog_hypo_chromHMM.txt	
grep 'Hypo' ../compare_marks/TE_combine_marks.txt > test
python ~/bin/TE_landscape/get_TE_sample.py test human_mouse_ortholog_hypo.txt human_mouse_ortholog_hypo_chromHMM.txt 8

# mm10 TEs also hypomethylated in human	
#TE_landscape/Mouse/mouse_human_ortholog_hypo_mm10.txt	
write.table(as.matrix(unique(human_mouse_orthologs_mm10_hypo[,9:15])),file="Mouse/mouse_human_ortholog_hypo_mm10.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

# Lifted over to mm9	
#TE_landscape/Mouse/mouse_human_ortholog_hypo_mm9.txt	
liftOver mouse_human_ortholog_hypo_mm10.txt /bar/genomes/mm10/chainFiles/mm10ToMm9.over.chain mouse_human_ortholog_hypo_mm9.txt out.txt #All successful

# Combined mm9 and mm10 coordinates	 
#TE_landscape/Mouse/mouse_human_ortholog_hypo.txt	
cut -f1-3 mouse_human_ortholog_hypo_mm9.txt | paste - mouse_human_ortholog_hypo_mm10.txt > mouse_human_ortholog_hypo.txt

# Intersected with chromHMM (mm9)	
#TE_landscape/Mouse/mouse_human_ortholog_hypo_chromHMM.txt	
while read line; do bedtools intersect -wo -a mouse_human_ortholog_hypo.txt -b ../raw_data/mouse/chromHMM/$line\.bed | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> mouse_human_ortholog_hypo_chromHMM.txt; done < ../sample_lists/mouse_chromHMM_samples.txt

# chromHMM state for mm10 TEs hypo in humans (mm9)	
#TE_landscape/Mouse/mouse_human_ortholog_hypo_chromHMM_sum.txt	
awk -v OFS='\t' '{a[$4, $5, $6, $7, $8, $9, $10, $14, $21]+=$20;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], sep[8], sep[9], a[i];}}' mouse_human_ortholog_hypo_chromHMM.txt > mouse_human_ortholog_hypo_chromHMM_sum.txt

# Updates 10/3/2017
## Adding chromHMM to WGBS. Getting hg19 chromHMM. 
 awk -v OFS='\t' '{if($12 < 0.3) print $0}' ../compare_marks/TE_combine_marks.txt > TE_combine_marks_meth.txt #Only TEs hypomethylated
 python ~/bin/TE_landscape/get_TE_sample.py TE_combine_marks_meth.txt hg19_ortholog_hypo.txt hg19_ortholog_hypo_chromHMM.txt 8
## Getting mm9 chromHMM.
 liftOver mm10_ortholog_hypo.txt /bar/genomes/mm10/chainFiles/mm10ToMm9.over.chain mm10_to_mm9_ortholog_hypo.txt out.txt
 cut -f1-3 mm10_to_mm9_ortholog_hypo.txt | paste - mm10_ortholog_hypo.txt > mm9_ortholog_hypo.txt
 while read line; do bedtools intersect -wo -a mm9_ortholog_hypo.txt -b ../raw_data/mouse/chromHMM/$line\.bed | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> mm9_ortholog_hypo_chromHMM.txt; done < ../sample_lists/mouse_chromHMM_samples.txt
 awk -v OFS='\t' '{a[$4, $5, $6, $7, $8, $9, $10, $14, $21]+=$20;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], sep[8], sep[9], a[i];}}' mm9_ortholog_hypo_chromHMM.txt > mm10_ortholog_hypo_chromHMM.txt
