# Human-mouse hypomethylated orthologs, chromHMM
# 5/29/2017, 5/30/2017, 7/26/2017, 7/31/2017, 10/3/2017, 10/10/2017

# hg19 TEs hypomethylated in human and mouse	
#Mouse/hg19_ortholog_hypo.txt	

# mm10 TEs hypomethylated in human and mouse
#Mouse/mm10_ortholog_hypo.txt

# hg19 chromHMM
## Pull out human chromHMM state for hypomethylated hg19 TEs
awk -v OFS='\t' '{if($11 == "Hypo") print $0}' compare_marks/TE_combine_marks.txt > Mouse/TE_combine_marks_meth.txt

## Human chromHMM state for hg19 TEs hypomethylated in mouse/human
python ~/bin/TE_landscape/get_TE_sample.py Mouse/TE_combine_marks_meth.txt Mouse/hg19_ortholog_hypo.txt Mouse/hg19_ortholog_hypo_chromHMM.txt 8

# mm10 chromHMM
## Liftover to mm9
ml kentUCSC
liftOver Mouse/mm10_ortholog_hypo.txt /bar/genomes/mm10/chainFiles/mm10ToMm9.over.chain Mouse/mm10_to_mm9_ortholog_hypo.txt Mouse/out.txt

# Combine mm9 and mm10 coordinates	 
cut -f1-3 Mouse/mm10_to_mm9_ortholog_hypo.txt | paste - Mouse/mm10_ortholog_hypo.txt > Mouse/mm9_ortholog_hypo.txt

# Intersect with chromHMM (mm9)	
while read line; do bedtools intersect -wo -a Mouse/mm9_ortholog_hypo.txt -b raw_data/mouse/chromHMM/$line\.bed | awk -v OFS='\t' -v sample=$line '{print $0, sample}' - >> Mouse/mm9_ortholog_hypo_chromHMM.txt; done < sample_lists/mouse_chromHMM_samples.txt

# chromHMM state for mm10 TEs hypomethylated in human/mouse (mm9)	
awk -v OFS='\t' '{a[$4, $5, $6, $7, $8, $9, $10, $14, $21]+=$20;}END{for(i in a) {split (i, sep, SUBSEP); print sep[1], sep[2], sep[3], sep[4], sep[5], sep[6], sep[7], sep[8], sep[9], a[i];}}' Mouse/mm9_ortholog_hypo_chromHMM.txt > Mouse/mm10_ortholog_hypo_chromHMM.txt
