# Load shuffled TE potential matrices, DNase

library(reshape2)

# Number of overlaps between shuffled TEs, DNase peaks
shuffled_DNase_potential = lapply(list.files(path="DNase/shuffled/",pattern="_DNase_peaks.txt",full.names = TRUE),function(x) read.table(x,sep='\t'))
shuffled_DNase_potential = lapply(shuffled_DNase_potential, setNames, nm =c("chromosome","start","stop","subfamily","class","family","strand","Sample","Peaks","Overlap"))
print("Loaded DNase matrices")
shuffled_DNase_potential = lapply(shuffled_DNase_potential,function(x) dcast(x,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks"))
print("Reshaping matrices")
shuffled_DNase_potential = lapply(shuffled_DNase_potential,function(x) replace(x,is.na(x),0))                                  
print("Replacing NA with 0")

# Number of samples a TE overlaps a DNase peak
shuffled_DNase_potential = lapply(shuffled_DNase_potential,function(y) transform(y,Samples = apply(y,1,function(x) sum(as.numeric(x[8:60]) > 0))))
print("Calculating samples")

save(shuffled_DNase_potential,file="R_datasets/shuffled_DNase.RData")