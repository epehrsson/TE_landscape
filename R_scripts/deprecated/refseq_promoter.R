# Load matrix of unique Refseq promoters

promoters = read.table("/bar/epehrsson/genic_features/RefSeq/refseq_promoters_unique_std.txt",sep='\t')
colnames(promoters) = c("chr","start","stop","strand")
promoters$Length = promoters$stop-promoters$start