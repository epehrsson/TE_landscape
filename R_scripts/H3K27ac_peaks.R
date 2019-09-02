# Loads a data frame of TE x sample H3K27ac peak overlap, for all 98 samples with H3K27ac data ("TE_H3K27ac_peaks")

# Load data frame listing TE, sample, and number of overlapping peaks
# Restricted to peaks where the summit overlaps the TE
print("Load H3K27ac")

TE_H3K27ac_peaks = read.table(file="H3K27ac/true_summit/rmsk_TEother_H3K27ac_summit.txt",sep='\t',
                              col.names=c("chromosome","start","stop","subfamily","class","family","strand","Sample","Peaks"))

# Reformat data frame (rows: TEs, columns: samples, values: number of overlapping peaks)
# Replacing missing values with 0
print("Reformat H3K27ac")

TE_H3K27ac_peaks = dcast(TE_H3K27ac_peaks,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks")
TE_H3K27ac_peaks[is.na(TE_H3K27ac_peaks)] = 0

# Update class names
TE_H3K27ac_peaks$class_update = convert_class(TE_H3K27ac_peaks$class)

# Count the number of samples each TE overlaps a H3K27ac peak summit
print("Number of samples overlapping H3K27ac peak")

TE_H3K27ac_peaks$Samples = apply(TE_H3K27ac_peaks,1,function(x) sum(as.numeric(x[8:105]) > 0))

# Count the number of samples each TE overlaps a H3K27ac peak summit, no cancer cell lines
TE_H3K27ac_peaks$Samples_noCancer = apply(TE_H3K27ac_peaks[,8:105],1,function(x) sum(as.numeric(x[which(metadata[match(colnames(TE_H3K27ac_peaks)[8:105],metadata$Sample),]$Exclude == "Include")]) > 0))

# Save the matrix
save(TE_H3K27ac_peaks,file="R_datasets/TE_H3K27ac_peaks.RData")