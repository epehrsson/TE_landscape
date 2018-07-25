# Peak/state widths
# 4/27/2016, 7/31/2017, 8/3/2017

# ChromHMM, block size by state and sample
#TE_landscape/chromHMM/chromHMM_blocks.txt
while read line; do awk -v OFS='\t' -v sample=$line '{print $3-$2, $4, sample}' raw_data/chromHMM/$line\_15_coreMarks_mnemonics.bed >> chromHMM/chromHMM_blocks.txt; done < sample_lists/mnemonics.txt

## Number of blocks per sample
while read line; do awk -v sample=$line -v OFS='\t' '{state[$4]+=1}END{for(i in state){print sample, i, state[i]}}' raw_data/chromHMM/$line\_15_coreMarks_mnemonics.bed >> chromHMM/chromHMM_peaks.txt; done < sample_lists/mnemonics.txt

#TE_landscape/DNase/peak_widths.txt
while read line; do awk -v OFS='\t' -v sample=$line '{print $3-$2, sample}' ../raw_data/DNase/DNase_narrow_peaks/$line\-DNase.macs2.narrowPeak >> peak_widths.txt; done < ../sample_lists/DNase_samples.txt

#TE_landscape/raw_data/H3K27ac/H3K27ac_narrow_peaks/peak_widths.txt
while read line; do awk -v OFS='\t' -v sample=$line '{print $3-$2, sample}' $line\-H3K27ac.narrowPeak >> peak_widths.txt; done < ../../../sample_lists/H3K27ac_samples.txt
