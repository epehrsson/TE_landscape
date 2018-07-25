# TEs always in a state

# Find the nearest gene
sort -k1,1 -k2,2n always_TEs.bed > always_TEs_sorted.bed
sort -k1,1 -k2,2n ~/genic_features/RefSeq/refseq_genes.txt | bedtools closest -a always_TEs_sorted.bed -b - -D b -t all > always_TEs_genes.txt

# Find the distance to the nearest gene
awk -v FS='\t' '{if(($8 == "DNase") && ($21 > 5000)) print $0}' always_TEs_genes.txt | wc -l
