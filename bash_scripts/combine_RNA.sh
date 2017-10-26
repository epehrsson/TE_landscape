# Average RNA-seq expression
# 5/31/2017, 6/6/2017

# Average expression per exon, raw	
#/scratch/ecp/RNA_average/refseq_exons.txt.sorted_E*.*.wig_average.txt [113 files]	

# Average expression per exon, raw, combined samples	 
#TE_landscape/RNAseq/refseq_exons_average.txt	
cut -f1-6 refseq_exons.txt.sorted_E000.neg.wig_average.txt > refseq_exons_average.txt
while read line; do cut -f9 refseq_exons.txt.sorted_$line\_average.txt | paste refseq_exons_average.txt - > refseq_exons_average; mv refseq_exons_average refseq_exons_average.txt ; done < ../samples.txt
while read line; do cut -f9 refseq_exons.txt.sorted_$line\_average.txt | paste refseq_exons_average.txt - > refseq_exons_average; mv refseq_exons_average refseq_exons_average.txt ; done < ../samples_agnostic.txt

# Average expression per TE, raw	
#/scratch/ecp/RNA_average/rmsk_TEother.txt.sorted_E*.*.wig_average.txt [113 files]	

# Average expression per TE, raw, combined samples	 
#TE_landscape/RNAseq/rmsk_TEother_average.txt	
cut -f1-7 rmsk_TEother.txt.sorted_E000.neg.wig_average.txt > rmsk_TEother_average.txt
while read line; do cut -f10 rmsk_TEother.txt.sorted_$line\_average.txt | paste rmsk_TEother_average.txt - > rmsk_TEother_average; mv rmsk_TEother_average rmsk_TEother_average.txt ; done < ../samples.txt
while read line; do cut -f10 rmsk_TEother.txt.sorted_$line\_average.txt | paste rmsk_TEother_average.txt - > rmsk_TEother_average; mv rmsk_TEother_average rmsk_TEother_average.txt ; done < ../samples_agnostic.txt
