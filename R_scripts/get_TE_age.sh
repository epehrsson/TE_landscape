# Age by TE
# 11/3/2016, 11/4/2016, 2/2/2017

# JC evolutionary distance for rmsk repeats	 
#TE_landscape/age/rmsk.txt.gz.JCage	
python calcAge_generalized.py rmsk.txt.gz.JCage

# JC evolutionary matrix for TEs, standard chromosomes	 
#TE_landscape/age/rmsk_TE_JCage.txt	
awk -v OFS='\t' '{if(($6 == "LTR" || $6 == "DNA" || $6 == "SINE" || $6 == "LINE") && ($1 !~ /_/)) print $1, $2, $3, $4, $5, $6, $7, $8, $9}' rmsk.txt.gz.JCage > rmsk_TE_JCage.txt
#TE_landscape/age/rmsk_other_JCage.txt	
awk -v OFS='\t' '{if(($6 == "LTR?" || $6 == "DNA?" || $6 == "SINE?" || $6 == "LINE?" || $6 == "Unknown?" || $6 == "Unknown" || $6 == "Other" || $6 == "RC") && ($1 !~ /_/)) print $1, $2, $3, $4, $5, $6, $7, $8, $9}' rmsk.txt.gz.JCage > rmsk_other_JCage.txt
