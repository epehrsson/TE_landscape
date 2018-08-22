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
## Mouse chromHMM state for mm10 TEs hypomethylated in human/mouse
