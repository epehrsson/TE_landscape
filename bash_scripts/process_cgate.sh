# Get hg19 coordaintes for C-GATE TEs in hg18 coordinates
# 1/7/18
# C_gate_hg18 and C_gate_hg19.bed are pulled from the CGate Excel

 liftOver features/C_gate_hg18 /bar/genomes/hg18/chainFiles/hg18ToHg19.over.chain.gz features/C_gate_hg18.bed features/out.txt
 cat features/C_gate_hg19.bed features/C_gate_hg18.bed > features/C_gate.txt

