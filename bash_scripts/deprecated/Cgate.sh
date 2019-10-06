# All code related to C-Gate analysis

# CGate TE coordinates	
#TE_landscape/features/C_gate.txt	
https://sites.google.com/site/tecatalog/welcome/home

# Get hg19 coordaintes for C-GATE TEs in hg18 coordinates
# C_gate_hg18 and C_gate_hg19.bed are pulled from the CGate Excel
liftOver features/C_gate_hg18 /bar/genomes/hg18/chainFiles/hg18ToHg19.over.chain.gz features/C_gate_hg18.bed features/out.txt
cat features/C_gate_hg19.bed features/C_gate_hg18.bed > features/C_gate.txt

# CGate TE coordinates intersected with all TEs	 
#TE_landscape/features/intersect_features/CGate_TE.txt	
bedtools intersect -wo -a C_gate.txt -b rmsk_TE.txt > CGate_TE.txt
bedtools intersect -wo -a C_gate.txt -b rmsk_other.txt > CGate_other.txt

# Updated 1/7/18 after hg18->hg19 liftover
bedtools intersect -wo -a features/C_gate.txt -b features/TEs/rmsk_TEother.txt > features/intersect_features/CGate_TE.txt
