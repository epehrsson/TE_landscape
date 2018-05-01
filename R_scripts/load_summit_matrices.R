# Load matrices of TE x sample x DNase/H3K27ac peak overlap (summit only)

## DNase

TE_DNase_summit = read.table(file="DNase/rmsk_TEother_DNase_summit.txt",sep='\t')
colnames(TE_DNase_summit) = c("chromosome","start","stop","subfamily","class","family","strand","Sample","Peaks")

# Number of overlaps between TEs, DNase peaks
TE_DNase_summit = dcast(TE_DNase_summit,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks")
TE_DNase_summit[is.na(TE_DNase_summit)] = 0
TE_DNase_summit$class_update = convert_class(TE_DNase_summit$class)

# Number of samples a TE overlaps a DNase peak
TE_DNase_summit$Samples = apply(TE_DNase_summit,1,function(x) sum(as.numeric(x[8:60]) > 0))

# Number of samples a TE overlaps a DNase peak, no cancer cell lines
TE_DNase_summit$Samples_noCancer = apply(TE_DNase_summit[,8:60],1,function(x) sum(as.numeric(x[which(metadata[match(colnames(TE_DNase_summit)[8:60],metadata$Sample),]$Exclude == "Include")]) > 0))

## H3K27ac

TE_H3K27ac_summit = read.table(file="H3K27ac/rmsk_TEother_H3K27ac_summit.txt",sep='\t')
colnames(TE_H3K27ac_summit) = c("chromosome","start","stop","subfamily","class","family","strand","Sample","Peaks")

# Number of overlaps between TEs, H3K27ac peaks
TE_H3K27ac_summit = dcast(TE_H3K27ac_summit,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks")
TE_H3K27ac_summit[is.na(TE_H3K27ac_summit)] = 0
TE_H3K27ac_summit$class_update = convert_class(TE_H3K27ac_summit$class)

# Number of samples a TE overlaps a H3K27ac peak
TE_H3K27ac_summit$Samples = apply(TE_H3K27ac_summit,1,function(x) sum(as.numeric(x[8:105]) > 0))

# Number of samples a TE overlaps a H3K27ac peak, no cancer cell lines
TE_H3K27ac_summit$Samples_noCancer = apply(TE_H3K27ac_summit[,8:105],1,function(x) sum(as.numeric(x[which(metadata[match(colnames(TE_H3K27ac_summit)[8:105],metadata$Sample),]$Exclude == "Include")]) > 0))

save(list=c("TE_DNase_summit","TE_H3K27ac_summit"),file="R_datasets/TE_summit.RData")