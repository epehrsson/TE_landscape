# Load shuffled TE potential matrices, H3K27ac

# Number of overlaps between shuffled TEs, H3K27ac peaks
print("Load H3K27ac matrices")
shuffled_H3K27ac_potential = lapply(seq(1,10,1),function(x) read.table(paste("H3K27ac/shuffled/rmsk_TE_",x,"_H3K27ac_summit.txt",sep=""),sep='\t',
                                                                       col.names=c(TE_coordinates[c(1:4,6,5,7)],"Sample","Peaks")))

print("Reshape matrices")
shuffled_H3K27ac_potential = lapply(shuffled_H3K27ac_potential,function(x) dcast(x,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks"))

print("Replace NA with 0")
shuffled_H3K27ac_potential = lapply(shuffled_H3K27ac_potential,function(x) replace(x,is.na(x),0))                                  

# Number of samples a TE overlaps a H3K27ac peak
print("Calculate samples")
shuffled_H3K27ac_potential = lapply(shuffled_H3K27ac_potential,function(y) transform(y,Samples = apply(y,1,function(x) sum(as.numeric(x[8:105]) > 0))))

save(shuffled_H3K27ac_potential,file="R_datasets/shuffled_H3K27ac.RData")