# Load matrix of TE x sample x DNase peak overlap
# See 5/24/2017, 6/5/2017
# Updated with summit only 5/1/18

print("Load DNase")

TE_DNase_peaks = read.table(file="DNase/rmsk_TEother_DNase_summit.txt",sep='\t')
colnames(TE_DNase_peaks) = c("chromosome","start","stop","subfamily","class","family","strand","Sample","Peaks")

# Number of overlaps between TEs, DNase peaks
print("Reformat DNase")

TE_DNase_peaks = dcast(TE_DNase_peaks,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks")
TE_DNase_peaks[is.na(TE_DNase_peaks)] = 0
TE_DNase_peaks$class_update = convert_class(TE_DNase_peaks$class)

# Number of samples a TE overlaps a DNase peak
print("Number of samples overlapping DNase peak")

TE_DNase_peaks$Samples = apply(TE_DNase_peaks,1,function(x) sum(as.numeric(x[8:60]) > 0))

# Number of samples a TE overlaps a DNase peak, no cancer cell lines
TE_DNase_peaks$Samples_noCancer = apply(TE_DNase_peaks[,8:60],1,function(x) sum(as.numeric(x[which(metadata[match(colnames(TE_DNase_peaks)[8:60],metadata$Sample),]$Exclude == "Include")]) > 0))

save(TE_DNase_peaks,file="R_datasets/TE_DNase_peaks.RData")