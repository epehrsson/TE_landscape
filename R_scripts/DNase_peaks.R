# TE x sample x DNase peak overlap (form matrix)
# See 5/24/2017, 6/5/2017

# Overlaps between TEs, DNase peaks
TE_DNase_peaks = read.table(file="DNase_peaks/rmsk_TEother_DNase_peaks_filter.txt",sep='\t')
colnames(TE_DNase_peaks) = c("chromosome","start","stop","subfamily","class","family","strand","Sample","Peaks","Overlap")

# Overlap length between TEs and DNase peaks
TE_DNase_peaks_overlap = dcast(TE_DNase_peaks,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Overlap")
TE_DNase_peaks_overlap[is.na(TE_DNase_peaks_overlap)] = 0

# Number of overlaps between TEs and DNase peaks
TE_DNase_peaks = dcast(TE_DNase_peaks,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks")
TE_DNase_peaks[is.na(TE_DNase_peaks)] = 0

# Number of samples a TE overlaps a DNase peak
TE_DNase_peaks$Samples = apply(TE_DNase_peaks,1,function(x) sum(as.numeric(x[8:60]) > 0))
TE_DNase_peaks$class_update = TE_DNase_peaks$class
TE_DNase_peaks$class_update = factor(TE_DNase_peaks$class_update,levels=c("DNA","LINE","LTR","SINE","Other","RC","Unconfident"))
TE_DNase_peaks[which(TE_DNase_peaks$class %in% c("DNA?","LINE?","LTR?","SINE?","Unknown","Unknown?")),]$class_update = "Unconfident"

# Potential
mean(TE_DNase_peaks$Samples)
mean(c(TE_DNase_peaks$Samples,rep(0,2589804))) 