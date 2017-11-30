# Run HOMER on enrichment candidates
input_list=$1 #enrichment/candidate_enrich_NEW_subfamily.txt

# HOMER results by subfamily (in state when enriched vs. never in state)
while read subfamily state
do 
  echo $subfamily $state
  findMotifsGenome.pl enrichment/$subfamily\_$state\_enriched.bed hg19 enrichment/HOMER_$subfamily\_$state\_enriched -size given -nomotif -bg enrichment/$subfamily\_no$state\.bed
  mv enrichment/HOMER_$subfamily\_$state\_enriched/knownResults.txt enrichment/HOMER_$subfamily\_$state\_enriched_knownResults.txt
  rm -r enrichment/HOMER_$subfamily\_$state\_enriched
done < $input_list 2>> enrichment/HOMER_output_NEW.txt

# Combine results
for file in enrichment/HOMER_*_enriched_knownResults.txt
do 
  awk -v OFS='\t' -v a=$file '{print a, $0}' $file >> enrichment/HOMER_enriched_knownResults.txt
done

rm enrichment/HOMER_*_enriched_knownResults.txt
