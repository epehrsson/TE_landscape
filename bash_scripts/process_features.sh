# Create new features for overlap
# 5/3/2016, 6/27/2016, 8/26/2016, 1/18/2017, 1/19/2017, 5/25/2017, 5/30/2017, 6/15/2017, 8/7/2017, 8/23/2017, 8/25/2017

# Genome
# File of chromosome sizes from UCSC, filtered to chr 1-22, X, Y	 
#TE_landscape/features/hg19_standard.genome	
head -n25 hg19.genome > hg19_standard.genome

# Segwey promoters

# Unique promoters	 
#TE_landscape/features/Segway_promoters/promoters.txt	
awk -v OFS='\t' '{print $1, $2, $3}' hg19_promoters.txt | sort | uniq > promoters.txt
sed 's/^/chr/' #To add "chr" to beginning of each line
#TE_landscape/features/Segway_promoters/promoters_merge.txt		
sort -k1,1V -k2,2n -k3,3n promoters.txt > promoters_sort.txt
bedtools merge -i promoters_sort.txt

# Refseq

# Expand 500bp downstream of RefSeq gene TSS	 
#genic_features/refseq_promoters.txt	
epehrsson@chimay:~/genic_features$ bedtools slop -i refseq_up2000.txt -g hg19.genome -l 0 -r 500 -s > refseq_promoters.txt

# Intergenic regions (no gene)	 
#genic_features/refseq_intergenic.txt	
bedtools complement -i refseq_genes.txt -g hg19.genome.standard > refseq_intergenic.txt

# RefSeq features, coding	 
#genic_features/refseq_[feature]_pc.txt [7 files]	
for file in refseq*.txt; do awk -v OFS='\t' '{if(substr($4,0,2) == "NM") print $0}' $file > $(basename "$file" .txt)\_pc.txt; done            
# RefSeq features, non-coding	 
#genic_features/refseq_[feature]_nc.txt [6 files]	
for file in refseq*.txt; do awk -v OFS='\t' '{if(substr($4,0,2) == "NR") print $0}' $file > $(basename "$file" .txt)\_nc.txt; done

# Unique Refseq promoter coordinates	 
#genic_features/RefSeq/refseq_promoters_unique.txt	
awk -v OFS='\t' '{print $1, $2, $3, $6}' ~/genic_features/RefSeq/refseq_promoters.txt | sort | uniq > refseq_promoters_unique.txt
# Filtered to standard chromosomes	 
#genic_features/RefSeq/refseq_promoters_unique_std.txt	
awk -v OFS='\t' '{if($1 !~ /_/) print $0}' refseq_promoters_unique.txt > refseq_promoters_unique_std.txt

# Number overlapping each other
 sort -k1,1V -k2,2n -k3,3n ~/genic_features/RefSeq/refseq_promoters_unique_std.txt | bedtools merge -i - > test.txt

# Unique Refseq exon coordinates, filtered to standard chromosomes
 awk -v OFS='\t' '{if($1 !~ /_/) print $1, $2, $3, $6}' ~/genic_features/RefSeq/refseq_exons.txt | sort | uniq > ~/genic_features/RefSeq/refseq_exons_unique.txt

# Number overlapping each other
 sort -k1,1V -k2,2n -k3,3n ~/genic_features/RefSeq/refseq_exons_unique.txt | bedtools merge -i - > test.txt

# GENCODE

# Extended to 500bp downstream	 
#genic_features/Gencodev19_promoter.txt	
bedtools slop -i Gencodev19_up2000.txt -g hg19.genome -l 0 -r 500 -s > Gencodev19_promoter.txt

# 2000 bp upstream, 500bp downstream of gene body	 
#genic_features/Gencode_v19_genes_margin.txt	
bedtools slop -i Gencode_v19_genes.txt -g hg19.genome -l 2000 -r 500 -s > Gencode_v19_genes_margin.txt

# Shuffled

# hg19 gaps, sorted	 
#TE_landscape/features/shuffled_TEs/gap_sorted.txt	
tail -n+2 /bar/genomes/hg19/gap/gap.txt | awk -v OFS='\t' '{print $2, $3, $4}' - | sort -k1,1V -k2,2n -k3,3n - > gap_sorted.txt

# Shuffled TE positions (no gaps)	 
#TE_landscape/features/shuffled_TEs/rmsk_TE_shuffle_#.txt [10 files]	
for i in {1..10}; do bedtools shuffle -i ~/TE_landscape/features/TEs/rmsk_TEother.txt -g ~/TE_landscape/features/hg19_standard.genome -excl gap_sorted.txt > rmsk_TE_shuffle_$i\.txt; done
