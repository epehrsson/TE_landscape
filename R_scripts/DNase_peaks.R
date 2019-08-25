# Loads a data frame of TE x sample DHS peak overlap, for all 53 samples with DHS data ("TE_DNase_peaks")

# Load data frame listing TE, sample, and number of overlapping peaks
# Restricted to peaks where the summit overlaps the TE
print("Load DNase")

TE_DNase_peaks = read.table(file="DNase/true_summit/rmsk_TEother_DNase_summit.txt",sep='\t')
colnames(TE_DNase_peaks) = c("chromosome","start","stop","subfamily","class","family","strand","Sample","Peaks")

# Reformat data frame (rows: TEs, columns: samples, values: number of overlapping peaks)
# Replacing missing values with 0
print("Reformat DNase")

TE_DNase_peaks = dcast(TE_DNase_peaks,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks")
TE_DNase_peaks[is.na(TE_DNase_peaks)] = 0

# Update class names
TE_DNase_peaks$class_update = convert_class(TE_DNase_peaks$class)

# Count the number of samples each TE overlaps a DHS peak summit
print("Number of samples overlapping DNase peak")

TE_DNase_peaks$Samples = apply(TE_DNase_peaks,1,function(x) sum(as.numeric(x[8:60]) > 0))

# Count the number of samples each TE overlaps a DHS peak summit, no cancer cell lines
TE_DNase_peaks$Samples_noCancer = apply(TE_DNase_peaks[,8:60],1,function(x) sum(as.numeric(x[which(metadata[match(colnames(TE_DNase_peaks)[8:60],metadata$Sample),]$Exclude == "Include")]) > 0))

# Save the matrix
save(TE_DNase_peaks,file="R_datasets/TE_DNase_peaks.RData")