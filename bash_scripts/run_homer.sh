# Run HOMER on enrichment candidates
input_list=$1 # subfamily\tstate

# HOMER results by subfamily (in state when enriched vs. never in state)
while read subfamily state
do 
  echo $subfamily $state
  findMotifsGenome.pl enrichment/bedfiles/$subfamily\_$state\_enriched.bed hg19 enrichment/homer/HOMER_$subfamily\_$state\_enriched -size given -nomotif -bg enrichment/bedfiles/$subfamily\_never_$state\.bed
  awk -v OFS='\t' -v subfam=$subfamily -v state=$state '{print subfam, state, $0}' enrichment/homer/HOMER_$subfamily\_$state\_enriched/knownResults.txt >> enrichment/homer/HOMER_enriched_knownResults.txt
done < $input_list 2>> enrichment/homer/HOMER_output.txt
