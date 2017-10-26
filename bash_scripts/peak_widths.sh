# Peak/state widths
# 4/27/2016, 7/31/2017, 8/3/2017

#TE_landscape/chromHMM/chromHMM_blocks.txt
for file in chromHMM_bedfiles/*.bed; do awk '{print ($3-$2), $4}' $file >> chromHMM_blocks.txt ; done

#TE_landscape/DNase/peak_widths.txt
while read line; do awk -v OFS='\t' -v sample=$line '{print $3-$2, sample}' ../raw_data/DNase/DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak >> peak_widths.txt; done < ../sample_lists/DNase_samples.txt

#TE_landscape/raw_data/H3K27ac/H3K27ac_narrow_peaks/peak_widths.txt
while read line; do awk -v OFS='\t' -v sample=$line '{print $3-$2, sample}' $line\-H3K27ac.narrowPeak >> peak_widths.txt; done < ../../../sample_lists/H3K27ac_samples.txt
