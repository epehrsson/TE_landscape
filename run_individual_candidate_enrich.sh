# Combined stats for subfamilies enriched in a state in at least two samples
 while read state; do awk -v OFS='\t' -v state=$state '{print state, $0}' enrichment/candidate_$state\_stat.txt >> enrichment/candidate_stats.txt; done < enrichment/states.txt

# Bedfiles of subfamily members in state in samples where subfamily is enriched for all candidate lists
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_1_TssA_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_2_TssAFlnk_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_3_TxFlnk_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_4_Tx_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_6_EnhG_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_7_Enh_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_8_ZNF.Rpts_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_9_Het_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_10_TssBiv_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_11_BivFlnk_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_12_EnhBiv_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_13_ReprPC_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_DNase_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_H3K27ac_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_Hypomethylated_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_Intermediate_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_Hypermethylated_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_Missing_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_5_TxWk_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_14_ReprPCWk_enriched.txt
python ~/bin/TE_landscape/get_enriched_TEs.py enrichment/candidate_15_Quies_enriched.txt

while read state; do awk -v OFS='\t' -v state=$state '{print state, $0}' enrichment/*_$state\_enrich*.bed >> enrichment/state_enriched.bed; done < enrichment/states.txt
wc -l enrichment/*_*_enriched.bed > enrichment/subfamily_state_enriched_counts.txt #Number of members in state when enriched
