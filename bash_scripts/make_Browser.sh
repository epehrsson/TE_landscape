 for file in raw_data/mouse/chromHMM/*.bed; do sort -k1,1 -k2,2n $file > ~/public_html/mouseENCODE/$(basename $file); bgzip ~/public_html/mouseENCODE/$(basename $file); done
 for file in ~/TE_landscape/raw_data/mouse/WGBS/*.bed; do awk -v OFS='\t' '{print $1, $2, $3, $10}' $file | sort -k1,1 -k2,2n - > $(basename $file .bed)_readDepth.bed; bgzip ~/public_html/mouseENCODE/$(basename $file .bed)_readDepth.bed; done
 for file in ~/TE_landscape/raw_data/mouse/WGBS/*.bed; do awk -v OFS='\t' '{print $1, $2, $3, $11}' $file | sort -k1,1 -k2,2n - > $(basename $file .bed)_meth.bed; bgzip ~/public_html/mouseENCODE/$(basename $file .bed)_meth.bed; done
 for file in *.bed.gz; do tabix -p bed $file; done #In public_html/mouseENCODE

