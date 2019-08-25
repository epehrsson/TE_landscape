# Loads data frames of shuffled TE x sample H3K27ac peak overlap, for all 98 samples with H3K27ac data ("shuffled_H3K27ac_potential")
# For 10 iterations of shuffled TEs

# Load data frames listing shuffled TE, sample, and number of overlapping peaks
# Restricted to peaks where the summit overlaps the TE
print("Load H3K27ac matrices")
shuffled_H3K27ac_potential = lapply(seq(1,10,1),function(x) read.table(paste("H3K27ac/shuffled/true_summit/rmsk_TE_",x,"_H3K27ac_summit.txt",sep=""),sep='\t',
                                                                       col.names=c(TE_coordinates[c(1:4,6,5,7)],"Sample","Peaks")))

# Reformat data frame (rows: TEs, columns: samples, values: number of overlapping peaks)
print("Reshape matrices")
shuffled_H3K27ac_potential = lapply(shuffled_H3K27ac_potential,function(x) dcast(x,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks"))

# Replace missing values with 0
print("Replace NA with 0")
shuffled_H3K27ac_potential = lapply(shuffled_H3K27ac_potential,function(x) replace(x,is.na(x),0))                                  

# Count the number of samples each TE overlaps a H3K27ac peak summit
print("Calculate samples")
shuffled_H3K27ac_potential = lapply(shuffled_H3K27ac_potential,function(y) transform(y,Samples = apply(y,1,function(x) sum(as.numeric(x[8:105]) > 0))))

# Save the matrices
save(shuffled_H3K27ac_potential,file="R_datasets/shuffled_H3K27ac.RData")