# GENCODE overlap

state=$1

# Overlap of TEs in state when enriched with GENCODE promoters (2000bp upstream, 500bp downstream TSS)
for file in enrichment/*_$state\_enriched.bed; do sort -k1,1V -k2,2 $file | bedtools intersect -wo -a - -b ~/genic_features/Gencode/Gencodev19_promoter.txt >> enrichment/Gencodev19_promoter_$state\_enriched.txt; done

# Overlap of TEs in state when enriched with GENCODE genes (2000 bp upstream, 500bp downstream gene)
for file in enrichment/*_$state\_enriched.bed; do sort -k1,1V -k2,2 $file | bedtools intersect -wo -a - -b ~/genic_features/Gencode/Gencode_v19_genes_margin.txt >> enrichment/Gencodev19_genes_$state\_enriched.txt; done
