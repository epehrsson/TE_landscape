# Human-mouse orthologs, chromHMM
# 5/29/2017, 5/30/2017, 7/26/2017, 7/31/2017, 10/3/2017, 10/10/2017

# hg19 chromHMM
## chromHMM state for all orthologous hg19 TEs
## Confirmed that it only pulled out TEs with matching subfamily and that the number of unique TEs is 269096 for all
cut -f5-11 Mouse/liftover/hg19_mm10_TE_intersect_same.bed | uniq | sort -k1,1 -k2,2n > Mouse/hg19_orthologs.txt
for file in E066 E071 E081 E083 E084 E086 E088 E092 E094 E096 E105 E106; do tail -n +2 /scratch/ecp/pandas/$file | sort -k1,1 -k2,2n - | bedtools intersect -sorted -wo -a Mouse/hg19_orthologs.txt -b - -f 1 -r > Mouse/hg19_ortholog_chromHMM_$file\.txt; done
for file in Mouse/chromHMM/hg19_ortholog_chromHMM_*.txt; do cut -f8-18 $file >> Mouse/chromHMM/hg19_chromHMM_TE.txt; done

