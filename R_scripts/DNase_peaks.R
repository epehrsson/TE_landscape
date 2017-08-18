# Load matrix of TE x sample x DNase peak overlap
# See 5/24/2017, 6/5/2017

library(reshape2)

# Number of overlaps between TEs, DNase peaks
TE_DNase_peaks = read.table(file="DNase/rmsk_TEother_DNase_peaks_filter.txt",sep='\t')
colnames(TE_DNase_peaks) = c("chromosome","start","stop","subfamily","class","family","strand","Sample","Peaks","Overlap")
TE_DNase_peaks = dcast(TE_DNase_peaks,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks")
TE_DNase_peaks[is.na(TE_DNase_peaks)] = 0

TE_DNase_peaks$class_update = TE_DNase_peaks$class
TE_DNase_peaks$class_update = factor(TE_DNase_peaks$class_update,levels=c("DNA","LINE","LTR","SINE","SVA","Other"))
TE_DNase_peaks[which(TE_DNase_peaks$class == "Other"),]$class_update = "SVA"
TE_DNase_peaks[which(TE_DNase_peaks$class %in% c("DNA?","LINE?","LTR?","SINE?","Unknown","Unknown?","RC")),]$class_update = "Other"

# Number of samples a TE overlaps a DNase peak
TE_DNase_peaks$Samples = apply(TE_DNase_peaks,1,function(x) sum(as.numeric(x[8:60]) > 0))

save(TE_DNase_peaks,file="R_scripts/TE_DNase_peaks.RData")
