# Mouse DNase analysis
# See 6/12/2017

# Overlap of mm10 TEs, DNase peaks
mm10_DNase = read.table("/scratch/ecp/DNase_mm10/mm10_orthologs_DNase_sum.txt",sep='\t')
colnames(mm10_DNase) = c(colnames(human_mouse_orthologs_mm10)[1:7],"File","Overlap","Peaks")

# Number of DNase peaks per TE 
mm10_DNase_peaks = dcast(data=mm10_DNase,mouse_chr_mm10+mouse_start_mm10+mouse_stop_mm10+mouse_subfamily+mouse_class+mouse_family+mouse_strand_mm10~File,value.var = "Peaks")
mm10_DNase_peaks[is.na(mm10_DNase_peaks)] = 0
