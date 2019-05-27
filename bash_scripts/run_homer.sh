# Run HOMER on candidate subfamilies

# List of candidates, format: subfamily\tstate
input_list=$1

# Identifiy transcription factor binding motifs (known only) encoded by subfamily members
while read subfamily state
do 
  echo $subfamily $state
  # Run on members in state when subfamily is enriched vs. members never in state as background
  findMotifsGenome.pl enrichment/bedfiles/$subfamily\_$state\_enriched.bed hg19 enrichment/homer/HOMER_$subfamily\_$state\_enriched -size given -nomotif -bg enrichment/bedfiles/$subfamily\_never_$state\.bed

  # Combine all results into one file
  awk -v OFS='\t' -v subfam=$subfamily -v state=$state '{print subfam, state, $0}' enrichment/homer/HOMER_$subfamily\_$state\_enriched/knownResults.txt >> enrichment/homer/HOMER_enriched_knownResults.txt
done < $input_list 2>> enrichment/homer/HOMER_output.txt
