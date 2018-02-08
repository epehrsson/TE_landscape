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

# For TEs in state when the sample is enriched, samples ever in state (1/12/18-1/15/18)
 python ~/bin/TE_landscape/pull_individual_TEs_state.py enrichment/LFSINE_Vert_3_TxFlnk_enriched.bed chromHMM/all_chromHMM_TE_sorted.txt LFSINE_Vert_3_TxFlnk_enriched_chromHMM.txt 3_TxFlnk 8
 python ~/bin/TE_landscape/pull_individual_TEs_state.py enrichment/MER31-int_3_TxFlnk_enriched.bed chromHMM/all_chromHMM_TE_sorted.txt enrichment/MER31-int_3_TxFlnk_enriched_chromHMM.txt 3_TxFlnk 8
 python ~/bin/TE_landscape/pull_individual_TEs_state.py enrichment/PRIMA41-int_3_TxFlnk_enriched.bed chromHMM/all_chromHMM_TE_sorted.txt enrichment/PRIMA41-int_3_TxFlnk_enriched_chromHMM.txt 3_TxFlnk 8
 for subfamily in MER68B Charlie10b Tigger4a; do python ~/bin/TE_landscape/pull_individual_TEs_state.py enrichment/$subfamily\_4_Tx_enriched.bed chromHMM/all_chromHMM_TE_sorted.txt enrichment/$subfamily\_4_Tx_enriched_chromHMM.txt 4_Tx 8; done
 for subfamily in HERV15-int HERV3-int HERVKC4-int HERVS71-int LTR10B1 LTR13 LTR25 LTR2C LTR30 LTR43B MER41G MER50C MER51A MER51C MER61F; do python ~/bin/TE_landscape/pull_individual_TEs_state.py enrichment/$subfamily\_1_TssA_enriched.bed chromHMM/all_chromHMM_TE_sorted.txt enrichment/$subfamily\_1_TssA_enriched_chromHMM.txt 1_TssA 8; done &
 awk '{if(($5 == "Other") && ($8 == "5_TxWk")) print $0}' chromHMM/all_chromHMM_other_sorted.txt > enrichment/SVA_5_TxWk_enriched_chromHMM.txt
 python ~/bin/TE_landscape/pull_individual_TEs_state.py enrichment/Tigger1a_Mars_5_TxWk_enriched.bed chromHMM/all_chromHMM_TE_sorted.txt enrichment/Tigger1a_Mars_5_TxWk_enriched_chromHMM.txt 5_TxWk 8
 for subfamily in Charlie10a Charlie10b L1P3b LTR71B MER68B; do python ~/bin/TE_landscape/pull_individual_TEs_state.py enrichment/$subfamily\_6_EnhG_enriched.bed chromHMM/all_chromHMM_TE_sorted.txt enrichment/$subfamily\_6_EnhG_enriched_chromHMM.txt 6_EnhG 8; done
 for subfamily in MER129 MER130 UCON12 UCON4 UCON6; do python ~/bin/TE_landscape/pull_individual_TEs_state.py enrichment/$subfamily\_7_Enh_enriched.bed chromHMM/all_chromHMM_other_sorted.txt enrichment/$subfamily\_7_Enh_enriched_chromHMM.txt 7_Enh 8; done
 while read line; do python ~/bin/TE_landscape/pull_individual_TEs_state.py enrichment/$line\_7_Enh_enriched.bed chromHMM/all_chromHMM_TE_sorted.txt enrichment/$line\_7_Enh_enriched_chromHMM.txt 7_Enh 8; done < enhancers.txt
 for subfamily in AmnSINE1 MER57E3 MER51A HERV15-int LTR10B1 LTR2C LTR43B LTR35B LTR26E LTR46 LTR77 LTR44 LTR4 LTR35A LTR3B MER88 MamGyp-int MamRep488 MER94B Charlie11 UCON29; do python ~/bin/TE_landscape/pull_individual_TEs_state.py enrichment/$subfamily\_2_TssAFlnk_enriched.bed chromHMM/all_chromHMM_TE_sorted.txt enrichment/$subfamily\_2_TssAFlnk_enriched_chromHMM.txt 2_TssAFlnk 8; done
 python ~/bin/TE_landscape/pull_individual_TEs_state.py enrichment/UCON26_2_TssAFlnk_enriched.bed chromHMM/all_chromHMM_other_sorted.txt enrichment/UCON26_2_TssAFlnk_enriched_chromHMM.txt 2_TssAFlnk 8
#Note, these are average methylation per TE, not overlap with CpGs in state. 
 python ~/bin/TE_landscape/pull_individual_TEs_state.py enrichment/UCON31_Hypomethylated_enriched.bed WGBS/TE_WGBS_state_sorted.txt enrichment/UCON31_Hypomethylated_enriched_WGBS.txt Hypomethylated 10
 python ~/bin/TE_landscape/pull_individual_TEs_state.py enrichment/MER57E3_Hypomethylated_enriched.bed WGBS/TE_WGBS_state_sorted.txt enrichment/MER57E3_Hypomethylated_enriched_WGBS.txt Hypomethylated 10

