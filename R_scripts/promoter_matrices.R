# Loads data frames of promoter x number of samples in each state
# For chromHMM state, average methylation, and overlap with DHS and H3K27ac peaks
## chromHMM_promoter_state: promoters, chromHMM state
## promoter_meth_average: promoters, average methylation level/WGBS state
## promoter_DNase_peaks: promoters, DHS peak overlap
## promoter_H3K27ac_peaks: promoters, H3K27ac peak overlap

# Data frame of unique promoter locations, standard chromosomes only
refseq_promoters = read.table("~/genic_features/RefSeq/refseq_promoters_unique_std.txt",sep='\t',col.names=c("chromosome","start","stop","strand"))

# Data frame with number of samples each promoter (row) is in each chromHMM state (column)
chromHMM_promoter_state = read.table("chromHMM/potential/refseq_promoters_chromHMM_summit_potential.txt",sep='\t',header=TRUE)
colnames(chromHMM_promoter_state)[5:19] = chromHMM_states

# Data frame with average DNA methylation of each promoter (row) in each sample (column)
promoter_meth_average = read.table("WGBS/Refseq_promoters/refseq_promoter_unique_CpG_Meth_average.txt",sep='\t',header=TRUE)

## Number of samples each promoter is in each methylation state
promoter_meth_average$Hypomethylated = apply(promoter_meth_average,1,function(x) sum(as.numeric(x[5:41]) < 0.3,na.rm=TRUE))
promoter_meth_average$Hypermethylated = apply(promoter_meth_average,1,function(x) sum(as.numeric(x[5:41]) > 0.7,na.rm=TRUE))
promoter_meth_average$Intermediate = apply(promoter_meth_average,1,function(x) sum(as.numeric(x[5:41]) <= 0.7 & as.numeric(x[5:41]) >= 0.3,na.rm=TRUE))
promoter_meth_average$Missing = apply(promoter_meth_average,1,function(x) sum(is.na(x[5:41])))

# Data frame of promoter, sample, and number of overlapping DHS peaks
## Restricted to peaks where the summit overlaps the promoter
promoter_DNase_peaks = read.table(file="DNase/true_summit/promoters/refseq_promoter_unique_DNase_summit.txt",sep='\t')
colnames(promoter_DNase_peaks) = c("chromosome","start","stop","strand","Sample","Peaks")

## Reformat matrix (row: promoter, column: sample, value: number of overlapping peak summits)
promoter_DNase_peaks = dcast(promoter_DNase_peaks,chromosome+start+stop+strand~Sample,value.var="Peaks")

## Expand to inclue all promoters
promoter_DNase_peaks = merge(refseq_promoters,promoter_DNase_peaks,by=c("chromosome","start","stop","strand"),all.x=TRUE)
promoter_DNase_peaks[is.na(promoter_DNase_peaks)] = 0

## Number of samples each promoter overlaps a DHS peak summit
promoter_DNase_peaks$Samples = apply(promoter_DNase_peaks,1,function(x) sum(as.numeric(x[5:57]) > 0))

# Data frame of promoter, sample, and number of overlapping H3K27ac peaks
## Restricted to peaks where the summit overlaps the promoter
promoter_H3K27ac_peaks = read.table(file="H3K27ac/true_summit/promoters/refseq_promoter_unique_H3K27ac_summit.txt",sep='\t')
colnames(promoter_H3K27ac_peaks) = c("chromosome","start","stop","strand","Sample","Peaks")

## Reformat matrix (row: promoter, column: sample, value: number of overlapping peak summits)
promoter_H3K27ac_peaks = dcast(promoter_H3K27ac_peaks,chromosome+start+stop+strand~Sample,value.var="Peaks")

## Expand to inclue all promoters
promoter_H3K27ac_peaks = merge(refseq_promoters,promoter_H3K27ac_peaks,by=c("chromosome","start","stop","strand"),all.x=TRUE)
promoter_H3K27ac_peaks[is.na(promoter_H3K27ac_peaks)] = 0

## Number of samples each promoter overlaps a H3K27ac peak summit
promoter_H3K27ac_peaks$Samples = apply(promoter_H3K27ac_peaks,1,function(x) sum(as.numeric(x[5:102]) > 0))

# Save data frames
save(list=c("chromHMM_promoter_state","promoter_meth_average","promoter_DNase_peaks","promoter_H3K27ac_peaks"),file="R_datasets/promoter_matrices.RData")