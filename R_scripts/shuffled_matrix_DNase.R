# Load shuffled TE potential matrices, DNase

# Number of overlaps between shuffled TEs, DNase peaks
print("Load DNase matrices")
shuffled_DNase_potential = lapply(seq(1,10,1),function(x) read.table(paste("DNase/shuffled/rmsk_TE_",x,"_DNase_summit.txt",sep=""),sep='\t',
                                                                     col.names=c(TE_coordinates[c(1:4,6,5,7)],"Sample","Peaks")))

print("Reshape matrices")
shuffled_DNase_potential = lapply(shuffled_DNase_potential,function(x) dcast(x,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks"))

print("Replace NA with 0")
shuffled_DNase_potential = lapply(shuffled_DNase_potential,function(x) replace(x,is.na(x),0))                                  

# Number of samples a TE overlaps a DNase peak
print("Calculate samples")
shuffled_DNase_potential = lapply(shuffled_DNase_potential,function(y) transform(y,Samples = apply(y,1,function(x) sum(as.numeric(x[8:60]) > 0))))

save(shuffled_DNase_potential,file="R_datasets/shuffled_DNase.RData")