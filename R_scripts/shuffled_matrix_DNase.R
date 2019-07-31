# Loads data frames of shuffled TE x sample DHS peak overlap, for all 53 samples with DHS data
# For 10 iterations of shuffled TEs

# Load data frames listing shuffled TE, sample, and number of overlapping peaks
# Restricted to peaks where the summit overlaps the TE
print("Load DNase matrices")
shuffled_DNase_potential = lapply(seq(1,10,1),function(x) read.table(paste("DNase/shuffled/true_summit/rmsk_TE_",x,"_DNase_summit.txt",sep=""),sep='\t',
                                                                     col.names=c(TE_coordinates[c(1:4,6,5,7)],"Sample","Peaks")))

# Reformat data frame (rows: TEs, columns: samples, values: number of overlapping peaks)
print("Reshape matrices")
shuffled_DNase_potential = lapply(shuffled_DNase_potential,function(x) dcast(x,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks"))

# Replace missing values with 0
print("Replace NA with 0")
shuffled_DNase_potential = lapply(shuffled_DNase_potential,function(x) replace(x,is.na(x),0))                                  

# Count the number of samples each TE overlaps a DHS peak summit
print("Calculate samples")
shuffled_DNase_potential = lapply(shuffled_DNase_potential,function(y) transform(y,Samples = apply(y,1,function(x) sum(as.numeric(x[8:60]) > 0))))

# Save the matrices
save(shuffled_DNase_potential,file="R_datasets/shuffled_DNase.RData")