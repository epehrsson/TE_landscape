# Class and subfamily sizes
# 8/18/2017, 8/19/2017

# Class lengths with/without chrY	 
#TE_landscape/features/TEs/class/class_lengths.txt	
awk -v OFS='\t' '{a[$4]+=$3-$2}END{for(i in a){print i, a[i]}}' TEother_class_merge.txt > class_lengths.txt
#TE_landscape/features/TEs/class/class_lengths_noY.txt		 
awk -v OFS='\t' '{if($1 != "chrY") a[$4]+=$3-$2}END{for(i in a){print i, a[i]}}' TEother_class_merge.txt > class_lengths_noY.txt

# Subfamily lengths with/without chrY	 
#TE_landscape/features/TEs/subfamily/subfamily_lengths.txt	
awk -v OFS='\t' '{a[$4]+=$3-$2}END{for(i in a){print i, a[i]}}' TEother_subfamily_merge.txt > subfamily_lengths.txt
#TE_landscape/features/TEs/subfamily/subfamily_lengths_noY.txt		 
awk -v OFS='\t' '{if($1 != "chrY") a[$4]+=$3-$2}END{for(i in a){print i, a[i]}}' TEother_subfamily_merge.txt > subfamily_lengths_noY.txt
