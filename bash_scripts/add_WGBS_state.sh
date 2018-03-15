# TE WGBS state (for state switching)
# 7/21/2017, 7/27/2017, 8/24/2017, 8/27/2017

# Methylation level of each TE x sample	 
#TE_landscape/WGBS/TE_WGBS_state.txt	
# Removed header	

# Methylation level of each TE x sample, sorted	 
#TE_landscape/WGBS/TE_WGBS_state_sorted.txt	
sort -k1,1V -k2,2n -k3,3n -k4,4 -k8,8 WGBS/TE_WGBS_state.txt > WGBS/TE_WGBS_state_sorted.txt
# Adding methylation states	 
awk -v OFS='\t' '{if($9 == "NA") print $0, "Missing"; else if ($9 < 0.3) print $0, "Hypomethylated"; else if ($9 > 0.7) print $0, "Hypermethylated"; else if (($9 <= 0.7) && ($9 >= 0.3)) print $0, "Intermediate";}' TE_WGBS_state_sorted.txt > TE_WGBS_state_sorted

# Methylation level/state of each TE x sample, sorted, by class	 
#TE_landscape/WGBS/class/[class]_WGBS_state_sorted.txt [6 files]	
while read line ; do awk -v OFS='\t' -v class=$line '{if($5 == class) print $0}' TE_WGBS_state_sorted.txt > $line\_WGBS_state_sorted.txt; done < ../features/TEs/class/TE_class.txt
awk -v OFS='\t' '{if($5 == "Other") print $0}' TE_WGBS_state_sorted.txt > SVA_WGBS_state_sorted.txt
awk -v OFS='\t' '{if(($5 != "Other") && ($5 != "LINE") && ($5 != "SINE") && ($5 != "DNA") && ($5 != "LTR")) print $0}' TE_WGBS_state_sorted.txt > Other_WGBS_state_sorted.txt

# Methylation state by subfamily x state
 awk '{print > "WGBS/subfamily/by_state/"$4"_"$10".txt"}' WGBS/TE_WGBS_state_sorted.txt
