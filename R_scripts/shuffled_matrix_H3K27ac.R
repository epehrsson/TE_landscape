# Load shuffled TE potential matrices, H3K27ac

library(reshape2)

# Number of overlaps between shuffled TEs, H3K27ac peaks
shuffled_H3K27ac_potential = lapply(list.files(path="H3K27ac/shuffled/",pattern="_H3K27ac_peaks.txt",full.names = TRUE),function(x) read.table(x,sep='\t'))
shuffled_H3K27ac_potential = lapply(shuffled_H3K27ac_potential, setNames, nm =c("chromosome","start","stop","subfamily","class","family","strand","Sample","Peaks","Overlap"))
print("Loaded H3K27ac matrices")
shuffled_H3K27ac_potential = lapply(shuffled_H3K27ac_potential,function(x) dcast(x,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks"))
print("Reshaping matrices")
shuffled_H3K27ac_potential = lapply(shuffled_H3K27ac_potential,function(x) replace(x,is.na(x),0))                                  
print("Replacing NA with 0")

# Number of samples a TE overlaps a H3K27ac peak
shuffled_H3K27ac_potential = lapply(shuffled_H3K27ac_potential,function(y) transform(y,Samples = apply(y,1,function(x) sum(as.numeric(x[8:105]) > 0))))
print("Calculating samples")

save(shuffled_H3K27ac_potential,file="R_scripts/shuffled_H3K27ac.RData")