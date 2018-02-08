 while read state; do awk -v OFS='\t' -v state=$state '{print state, $0}' enrichment/*_$state\_enrich*.bed >> enrichment/state_enriched.bed; done < enrichment/states.txt 
 wc -l enrichment/*_*_enriched.bed > enrichment/subfamily_state_enriched_counts.txt #Number of members in state when enriched

