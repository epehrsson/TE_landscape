# Merged class and subfamily files
# 2/2/2017, 2/3/2017, 3/2/2017, 6/5/2017, 8/18/2017

# Merged file for each class	
#TE_landscape/features/TEs/class/rmsk_[class].txt [13 files]	
while read line; do awk -v OFS='\t' -v class=$line '{if($5 == class)print $0}' rmsk_other.txt | bedtools merge -i - > TE_classes/rmsk_$line\.txt; done < other_class.txt &
while read line; do awk -v OFS='\t' -v class=$line '{if($5 == class)print $0}' rmsk_TE.txt | bedtools merge -i - > TE_classes/rmsk_$line\.txt; done < TE_class.txt &
awk -v OFS='\t' '{if(($5=="LTR?")||($5=="DNA?")||($5=="LINE?")||($5=="SINE?")||($5=="Unknown?")||($5=="Unknown")) print $0}' ../rmsk_other.txt | bedtools merge -i - > rmsk_Unconfident.txt
# Adding RC to Unconfident class	 
#TE_landscape/features/TEs/class/rmsk_Unconfident_RC.txt	
awk -v OFS='\t' '{if(($5=="LTR?")||($5=="DNA?")||($5=="LINE?")||($5=="SINE?")||($5=="Unknown?")||($5=="Unknown")||($5=="RC")) print $0}' ../rmsk_TEother.txt | bedtools merge -i - > rmsk_Unconfident_RC.txt
    
# Merged TE bases by class	 
#TE_landscape/features/TEs/class/TEother_class_merge.txt	
while read line; do awk -v OFS='\t' -v class=$line '{print $0, class}' rmsk_$line\.txt >> TEother_class_merge.txt ; done < TEother_class.txt
# Adding merged Unconfident-RC bases	 
awk -v OFS='\t' '{print $0, "Unconfident_RC"}' rmsk_Unconfident_RC.txt >> TEother_class_merge.txt
    
# Merged file for each subfamily	 
#TE_landscape/features/TEs/subfamily/TEother_subfamily_merge.txt	
while read line; do awk -v OFS='\t' -v subfam=$line '{if($4 == subfam)print $0}' ../rmsk_other.txt | bedtools merge -i - > rmsk_$line\.txt; done < other_subfamilies.txt
while read line; do awk -v OFS='\t' -v subfam=$line '{if($4 == subfam)print $0}' ../rmsk_TE.txt | bedtools merge -i - > rmsk_$line\.txt; done < subfamilies.txt
while read line; do awk -v OFS='\t' -v subfam=$line '{print $0, subfam}' rmsk_$line\.txt >> TEother_subfamily_merge.txt; done < other_subfamilies.txt
while read line; do rm rmsk_$line\.txt; done < other_subfamilies.txt
while read line; do awk -v OFS='\t' -v subfam=$line '{print $0, subfam}' rmsk_$line\.txt >> TEother_subfamily_merge.txt; rm rmsk_$line\.txt; done < subfamilies.txt
