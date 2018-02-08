bash bash_scripts/overlap_individual_TE_GENCODE.sh 1_TssA
bash bash_scripts/overlap_individual_TE_GENCODE.sh 2_TssAFlnk
bash bash_scripts/overlap_individual_TE_GENCODE.sh 3_TxFlnk
bash bash_scripts/overlap_individual_TE_GENCODE.sh 4_Tx
bash bash_scripts/overlap_individual_TE_GENCODE.sh 5_TxWk
bash bash_scripts/overlap_individual_TE_GENCODE.sh 6_EnhG
bash bash_scripts/overlap_individual_TE_GENCODE.sh 7_Enh
bash bash_scripts/overlap_individual_TE_GENCODE.sh Hypomethylated
bash bash_scripts/overlap_individual_TE_GENCODE.sh Intermediate
bash bash_scripts/overlap_individual_TE_GENCODE.sh DNase
bash bash_scripts/overlap_individual_TE_GENCODE.sh H3K27ac

# Combine results
while read state; do awk -v OFS='\t' -v state=$state '{print state, $0}' enrichment/Gencodev19_genes_$state\_enriched.txt >> enrichment/Gencodev19_genes_enriched.txt; done < enrichment/states.txt
while read state; do awk -v OFS='\t' -v state=$state '{print state, $0}' enrichment/Gencodev19_promoter_$state\_enriched.txt >> enrichment/Gencodev19_promoter_enriched.txt; done < enrichment/states.txt
