# Generates dataframes with the length of each chromosome annotated with each epigenetic state across samples,
# the proportion of the genome and individual TEs represented by each chromosome, 
# and the proportion of TEs on each chromosome that overlap each genic feature

## chrom_state_total - Length of each chromosome annotated with each chromHMM state, across all samples
## chr_state_WGBS_total - Number of CpGs annotated with each methylation state by chromosome, across all samples
## chr_state_DNase_total - Number of DHS peaks on each chromosome, across all samples
## chr_state_H3K27ac_total - Number of H3K27ac peaks on each chromosome, across all samples
## hg19_TEs - Length of each chromosome and the proportion of the genome represented by each chromosome, 
## plus the number/proportion of TEs on each chromosome
## hg19_features - Number/proportion of TEs on each chromosome overlapping each genic feature
## chromosome_potential - Z-score for each chromosome based on mean number of samples each TE is annotated with the state by chromosome


# Epigenetic state annotation by chromosome

## chromHMM: length of each chromosome annotated with each chromHMM state, by sample
chrom_state = read.table("chromHMM/chromosome_states.txt",sep='\t',col.names=c("Sample","chromosome","State","Bases"))
chrom_state$chromosome = factor(chrom_state$chromosome,levels=standard_chromosomes)
chrom_state$State = factor(chrom_state$State,levels=chromHMM_states)

## Length of each chromosome annotated with each chromHMM state, across all samples
chrom_state_total = ddply(chrom_state,.(chromosome,State),summarise,Bases=sum(as.numeric(Bases)))

## WGBS: number of CpGs annotated with each methylation state, by chromosome and sample
chr_state_WGBS = read.table("WGBS/chr_CpG_meth_states.txt",sep='\t',col.names=c("chromosome","Sample",meth_states))
chr_state_WGBS[is.na(chr_state_WGBS)] = 0
chr_state_WGBS$Sample = mapvalues(chr_state_WGBS$Sample,seq(4,40,1),as.vector(metadata[which(!is.na(metadata$WGBS)),]$Sample))
chr_state_WGBS$chromosome = factor(chr_state_WGBS$chromosome,levels=standard_chromosomes)
chr_state_WGBS = melt(chr_state_WGBS,id.vars=c("chromosome","Sample"),variable.name="State",value.name="CpGs")

## Number of CpGs annotated with each methylation state by chromosome, across all samples
chr_state_WGBS_total = ddply(chr_state_WGBS,.(chromosome,State),summarise,CpGs=sum(CpGs))

## DHS: Number of DHS peaks on each chromosome, by sample
chr_state_DNase = read.table("DNase/chr_DNase.txt",sep=" ",col.names=c("Sample","chromosome","Peaks"))
chr_state_DNase$chromosome = factor(chr_state_DNase$chromosome,levels=standard_chromosomes)

## Number of DHS peaks on each chromosome, across all samples
## Plus the length of the chromosome (bp)
chr_state_DNase_total = ddply(chr_state_DNase,.(chromosome),summarise,Peaks=sum(Peaks))
chr_state_DNase_total = merge(chr_state_DNase_total,hg19_genome,by.x="chromosome",by.y="chrom")

## H3K27ac: Number of H3K27ac peaks on each chromosome, by sample
chr_state_H3K27ac = read.table("H3K27ac/chr_H3K27ac.txt",sep=" ",col.names=c("Sample","chromosome","Peaks"))
chr_state_H3K27ac$chromosome = factor(chr_state_H3K27ac$chromosome,levels=standard_chromosomes)

## Number of H3K27ac peaks on each chromosome, across all samples
## Plus the length of the chromosome (bp)
chr_state_H3K27ac_total = ddply(chr_state_H3K27ac,.(chromosome),summarise,Peaks=sum(Peaks))
chr_state_H3K27ac_total = merge(chr_state_H3K27ac_total,hg19_genome,by.x="chromosome",by.y="chrom")


# Length of each chromosome (bp) and the proportion of the genome represented by each chromosome
# chr1-22, chrX, chrY, and chrM only
hg19_TEs = hg19_genome[which(hg19_genome$chrom %in% standard_chromosomes),]
hg19_TEs = ddply(hg19_TEs,.(),transform,Prop=size/sum(as.numeric(size)))[2:4]
hg19_TEs$chrom = factor(hg19_TEs$chrom,levels=standard_chromosomes)

# Add the number/proportion of TEs on each chromosome
hg19_TEs = merge(hg19_TEs,ddply(rmsk_TE,.(chromosome),summarize,TEs=length(chromosome)),by.x="chrom",by.y="chromosome",all.x=TRUE)
hg19_TEs = ddply(hg19_TEs,.(),transform,Prop_TEs=TEs/NUM_TE)[2:6]


# Number/proportion of TEs on each chromosome overlapping each genic feature
hg19_features = ddply(melt(rmsk_TE[,c("chromosome",cohorts)],id.var="chromosome",variable.name="Cohort",value.name="Overlap"),
             .(chromosome,Cohort),summarise,TEs_feature=sum(!is.na(Overlap)))
hg19_features = split_coding(hg19_features)
hg19_features = merge(hg19_TEs,hg19_features,by.x="chrom",by.y="chromosome")
hg19_features$Prop_TEs_feature = hg19_features$TEs_feature/hg19_features$TEs


# Mean number of samples each TE is annotated with each state, by chromosome
chromosome_potential = ddply(melt(rmsk_TE_measure[,c("chromosome",states)],id.var="chromosome",variable.name="State"),
                             .(chromosome,State),summarise,Samples=mean(na.omit(value)))
chromosome_potential$chromosome = factor(chromosome_potential$chromosome,levels=standard_chromosomes)
## For each state, calculate a Z-score for each chromosome based on mean number of samples each TE is annotated with the state
chromosome_potential = ddply(chromosome_potential,.(State),transform,Zscore=(Samples-mean(Samples))/sd(Samples))

rm(list=c("chrom_state", "chr_state_WGBS", "chr_state_DNase", "chr_state_H3K27ac"))