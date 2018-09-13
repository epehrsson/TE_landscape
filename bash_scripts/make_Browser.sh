# Process tracks for display on the Browser

# mm10 chromHMM categorical tracks
## Fourth column must be positive integer
while read a b c; do gunzip raw_data/mouse/chromHMM/$c\.bed.gz; awk -v OFS='\t' '{if($4 == "TssA"){print $1, $2, $3, 1} else if ($4 == "TssAFlnk1"){print $1, $2, $3, 2} else if ($4 == "TssAFlnk2"){print $1, $2, $3, 3} else if ($4 == "Tx1"){print $1, $2, $3, 4} else if ($4 == "Tx2"){print $1, $2, $3, 5} else if ($4 == "Enh"){print $1, $2, $3, 6} else if ($4 == "EnhLo1"){print $1, $2, $3, 7} else if ($4 == "EnhLo2"){print $1, $2, $3, 8} else if ($4 == "HetCons"){print $1, $2, $3, 9} else if ($4 == "HetFac"){print $1, $2, $3, 10} else if ($4 == "TssBiv"){print $1, $2, $3, 11} else if ($4 == "EnhPois1"){print $1, $2, $3, 12} else if ($4 == "EnhPois2"){print $1, $2, $3, 13} else if ($4 == "QuiesG"){print $1, $2, $3, 14} else if ($4 == "Quies"){print $1, $2, $3, 15}}' raw_data/mouse/chromHMM/$c\.bed | sort -k1,1 -k2,2n - > /tavern/epehrsson/mouseENCODE/signal_tracks/$c\.bed; bgzip /tavern/epehrsson/mouseENCODE/signal_tracks/$c\.bed; tabix -p bed /tavern/epehrsson/mouseENCODE/signal_tracks/$c\.bed.gz; done < Mouse/human_mouse_samples.txt &

# mm10 WGBS
for file in ~/TE_landscape/raw_data/mouse/WGBS/*.bed; do awk -v OFS='\t' '{print $1, $2, $3, $10}' $file | sort -k1,1 -k2,2n - > $(basename $file .bed)_readDepth.bed; bgzip ~/public_html/mouseENCODE/$(basename $file .bed)_readDepth.bed; done
for file in ~/TE_landscape/raw_data/mouse/WGBS/*.bed; do awk -v OFS='\t' '{print $1, $2, $3, $11}' $file | sort -k1,1 -k2,2n - > $(basename $file .bed)_meth.bed; bgzip ~/public_html/mouseENCODE/$(basename $file .bed)_meth.bed; done
for file in *.bed.gz; do tabix -p bed $file; done #In public_html/mouseENCODE

