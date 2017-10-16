# Process enrichment candidates
# 11/4/2016, 11/21/2016, 2/7/17, 2/25/2017, 5/23/2017, 5/24/2017, 5/25/2017, 5/30/2017, 5/31/2017

# GENCODE overlap

# Intersection of subfamily and RefSeq promoters	 
bedtools intersect -a LTR12E.txt -b ../../genic_features/refseq_promoters_merge.txt

# Intersection of 1_TssA candidate subfamilies and RefSeq promoters	 
while read line; do  bedtools intersect -wo -a enrichment/$line\.txt -b ~/genic_features/refseq_promoters_merge.txt | awk '{print $1, $2, $3}' - | sort | uniq | wc -l ; echo $line; done < candidate_1TssA.txt
while read line; do echo $line; tail -n +2 $line\_1TssA.bed | bedtools intersect -wo -a - -b ~/genic_features/refseq_promoters_merge.txt | awk '{print $1, $2, $3}' - | sort | uniq | wc -l ; done < candidate_1TssA.txt

# Overlap of TEs ever in 1_TssA state with GENCODE promoters	 
#TE_landscape/enrichment/Gencodev19_promoter_1TssA.txt	
for file in *_1TssA.bed; do tail -n +2 $file | bedtools intersect -wo -a - -b ~/genic_features/Gencodev19_promoter.txt >> Gencodev19_promoter_1TssA.txt; done

# Number of TEs with overlap 	 
awk '{print $1, $2, $3, $4}' Gencodev19_promoter_1TssA.txt | sort | uniq | awk '{a[$4]+=1}END{for(i in a){print i, a[i];}}' -

# Overlap of TEs in 1_TssA state when enriched with GENCODE promoters	 
#TE_landscape/enrichment/Gencodev19_promoter_1TssA_enrich.txt	
for file in *_enrich1TssA.bed; do sort -k1,1V -k2,2 $file | bedtools intersect -wo -a - -b ~/genic_features/Gencodev19_promoter.txt >> Gencodev19_promoter_1TssA_enrich.txt; done

# Overlap of TEs ever in state with GENCODE genes	 
#TE_landscape/enrichment/Gencodev19_genes_2TssAFlnk.txt	
for file in *_2TssAFlnk.bed; do bedtools intersect -wo -a $file -b ~/genic_features/Gencode_v19_genes.txt >> Gencodev19_genes_2TssAFlnk.txt; done
#TE_landscape/enrichment/Gencodev19_genes_3TxFlnk.txt		 
for file in *_3TxFlnk.bed; do bedtools intersect -wo -a $file -b ~/genic_features/Gencode_v19_genes.txt >> Gencodev19_genes_3TxFlnk.txt; done
#TE_landscape/enrichment/Gencodev19_genes_6EnhG.txt		 
for file in *_6EnhG.bed; do bedtools intersect -wo -a $file -b ~/genic_features/Gencode_v19_genes.txt >> Gencodev19_genes_6EnhG.txt; done

# Overlap of TEs in state when enriched with GENCODE genes (margin)	 
#TE_landscape/enrichment/Gencodev19_genes_2TssAFlnk_enrich.txt	
for file in *_enrich2TssAFlnk.bed; do sort -k1,1V -k2,2 $file | bedtools intersect -wo -a - -b ~/genic_features/Gencode_v19_genes_margin.txt >> Gencodev19_genes_2TssAFlnk_enrich.txt; done
#TE_landscape/enrichment/Gencodev19_genes_3TxFlnk_enrich.txt		 
for file in *_enrich3TxFlnk.bed; do sort -k1,1V -k2,2 $file | bedtools intersect -wo -a - -b ~/genic_features/Gencode_v19_genes_margin.txt >> Gencodev19_genes_3TxFlnk_enrich.txt; done
#TE_landscape/enrichment/Gencodev19_genes_6EnhG_enrich.txt		 
for file in *_enrich6EnhG.bed; do sort -k1,1V -k2,2 $file | bedtools intersect -wo -a - -b ~/genic_features/Gencode_v19_genes_margin.txt >> Gencodev19_genes_6EnhG_enrich.txt; done

# HOMER

# HOMER results by subfamily (all members)	 
#TE_landscape/enrichment/HOMER/HOMER_[subfamily]_knownResults.txt [157 files]	
while read line; do findMotifsGenome.pl $line\.txt hg19 HOMER_$line -size given -nomotif; done < enhancer_subfamilies.txt 2>> HOMER_enhancer_output.txt &

# HOMER output	
#TE_landscape/enrichment/HOMER/HOMER_enhancer_output.txt	

# Motifs enriched in MER57E3 copies in the 1_TssA state in at least 20 samples	 
#TE_landscape/enrichment/HOMER/HOMER_MER57E3_1TssA20_knownResults.txt	
findMotifsGenome.pl MER57E3_1TssA.bed hg19 HOMER_$line -size given -nomotif -bg MER57E3_no1TssA.bed

# HOMER results by subfamily (ever in state vs. never in state)	 
#TE_landscape/enrichment/HOMER/HOMER_*_1TssA_knownResults.txt [40 files]	
while read line; do findMotifsGenome.pl $line\_1TssA.bed hg19 HOMER_$line\_1TssA -size given -nomotif -bg $line\_no1TssA.bed; done < ../candidate_1TssA.txt
while read line; do findMotifsGenome.pl $line\_2TssAFlnk.bed hg19 HOMER_$line\_2TssAFlnk -size given -nomotif -bg $line\_no2TssAFlnk.bed; done < ../candidate_2TssAFlnk.txt &
while read line; do findMotifsGenome.pl $line\_3TxFlnk.bed hg19 HOMER_$line\_3TxFlnk -size given -nomotif -bg $line\_no3TxFlnk.bed; done < ../candidate_3TxFlnk.txt
while read line; do findMotifsGenome.pl $line\_6EnhG.bed hg19 HOMER_$line\_6EnhG -size given -nomotif -bg $line\_no6EnhG.bed; done < ../candidate_6EnhG.txt &

# HOMER results by subfamily (in state when enriched vs. never in state)	 
#TE_landscape/enrichment/HOMER/HOMER_*_enrich1TssA_knownResults.txt [39 files]	
while read line; do findMotifsGenome.pl $line\_enrich1TssA.bed hg19 HOMER_$line\_enrich1TssA -size given -nomotif -bg $line\_no1TssA.bed; done < candidate_1TssA.txt
#TE_landscape/enrichment/HOMER/HOMER_*_enrich2TssAFlnk_knownResults.txt [91 files]		 
while read line; do findMotifsGenome.pl $line\_enrich2TssAFlnk.bed hg19 HOMER_$line\_enrich2TssAFlnk -size given -nomotif -bg $line\_no2TssAFlnk.bed; done < ../candidate_2TssAFlnk.txt
#TE_landscape/enrichment/HOMER/HOMER_*_enrich3TxFlnk_knownResults.txt [13 files]		 
while read line; do findMotifsGenome.pl $line\_enrich3TxFlnk.bed hg19 HOMER_$line\_enrich3TxFlnk -size given -nomotif -bg $line\_no3TxFlnk.bed; done < ../candidate_3TxFlnk.txt
#TE_landscape/enrichment/HOMER/HOMER_*_enrich6EnhG_knownResults.txt [36 files]		 
while read line; do findMotifsGenome.pl $line\_enrich6EnhG.bed hg19 HOMER_$line\_enrich6EnhG -size given -nomotif -bg $line\_no6EnhG.bed; done < ../candidate_6EnhG.txt
#TE_landscape/enrichment/HOMER/HOMER_*_enrich7Enh_knownResults.txt [159 files]		 
while read line; do findMotifsGenome.pl $line\_enrich7Enh.bed hg19 HOMER_$line\_enrich7Enh -size given -nomotif -bg $line\_no7Enh.bed; done < ../candidate_7Enh.txt

# Combined 	 
#TE_landscape/enrichment/HOMER/HOMER_enrich_knownResults.txt	
while read line; do awk -v OFS='\t' -v a=$line '{print a, $0}' HOMER_$line\_knownResults.txt >> HOMER_enrich_knownResults.txt; done < enrich_results.txt
