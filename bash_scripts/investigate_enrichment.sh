# Process enrichment candidates
# 9/30/2016, 11/4/2016, 11/20/16, 11/21/2016, 2/7/17, 2/25/2017, 5/23/2017, 5/24/2017, 5/25/2017, 5/30/2017, 5/31/2017

# Individual subfamilies
# Subfamily coordinates in bed format	 
#TE_landscape/features/TEs/subfamily/individual/xx.txt [269 files]	
while read line; do grep -w $line ../rmsk_TE.txt | awk -v OFS='\t' '{print $1,$2,$3,$4,"0",$7}' - > $line\.txt; done < enhancer_subfamilies.txt

# Intersection of subfamily and RefSeq promoters	 
bedtools intersect -a LTR12E.txt -b ../../genic_features/refseq_promoters_merge.txt

# Intersection of 1_TssA candidate subfamilies and RefSeq promoters	 
while read line; do  bedtools intersect -wo -a enrichment/$line\.txt -b ~/genic_features/refseq_promoters_merge.txt | awk '{print $1, $2, $3}' - | sort | uniq | wc -l ; echo $line; done < candidate_1TssA.txt
while read line; do echo $line; tail -n +2 $line\_1TssA.bed | bedtools intersect -wo -a - -b ~/genic_features/refseq_promoters_merge.txt | awk '{print $1, $2, $3}' - | sort | uniq | wc -l ; done < candidate_1TssA.txt

# Number of TEs with overlap 	 
awk '{print $1, $2, $3, $4}' Gencodev19_promoter_1TssA.txt | sort | uniq | awk '{a[$4]+=1}END{for(i in a){print i, a[i];}}' -
