# Load matrix of TE x sample x H3K27ac peak overlap
# See 7/4/2017

library(reshape2)

# Number of overlaps between TEs, H3K27ac peaks
TE_H3K27ac_peaks = read.table(file="H3K27ac/rmsk_TEother_H3K27ac_peaks_filter.txt",sep='\t')
colnames(TE_H3K27ac_peaks) = c("chromosome","start","stop","subfamily","class","family","strand","Sample","Peaks","Overlap")
TE_H3K27ac_peaks = dcast(TE_H3K27ac_peaks,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks")
TE_H3K27ac_peaks[is.na(TE_H3K27ac_peaks)] = 0

TE_H3K27ac_peaks$class_update = TE_H3K27ac_peaks$class
TE_H3K27ac_peaks$class_update = factor(TE_H3K27ac_peaks$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
TE_H3K27ac_peaks[which(TE_H3K27ac_peaks$class == "Other"),]$class_update = "SVA"
TE_H3K27ac_peaks[which(TE_H3K27ac_peaks$class %in% c("DNA?","LINE?","LTR?","SINE?","Unknown","Unknown?","RC")),]$class_update = "Other"

# Number of samples a TE overlaps a H3K27ac peak
TE_H3K27ac_peaks$Samples = apply(TE_H3K27ac_peaks,1,function(x) sum(as.numeric(x[8:105]) > 0))

# Number of samples a TE overlaps a H3K27ac peak, no cancer cell lines
TE_H3K27ac_peaks$Samples_noCancer = apply(TE_H3K27ac_peaks[,8:105],1,function(x) sum(as.numeric(x[which(metadata[match(colnames(TE_H3K27ac_peaks)[8:105],metadata$Sample),]$Exclude == "Include")]) > 0))

save(TE_H3K27ac_peaks,file="R_datasets/TE_H3K27ac_peaks.RData")