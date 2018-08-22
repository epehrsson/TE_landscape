# Chromosome potential

# Overall epigenetic state by chromosome
chrom_state = read.table("chromHMM/chromosome_states.txt",sep='\t')
colnames(chrom_state) = c("Sample","chromosome","State","Bases")
chrom_state$chromosome = factor(chrom_state$chromosome,levels=standard_chromosomes)
chrom_state$State = factor(chrom_state$State,levels=chromHMM_states)
chrom_state_total = ddply(chrom_state,.(chromosome,State),summarise,Bases=sum(as.numeric(Bases)))

chr_state_WGBS = read.table("WGBS/chr_CpG_meth_states.txt",sep='\t')
colnames(chr_state_WGBS) = c("chromosome","Sample",meth_states)
chr_state_WGBS[is.na(chr_state_WGBS)] = 0
chr_state_WGBS$Sample = mapvalues(chr_state_WGBS$Sample,seq(4,40,1),as.vector(metadata[which(!is.na(metadata$WGBS)),]$Sample))
chr_state_WGBS$chromosome = factor(chr_state_WGBS$chromosome,levels=standard_chromosomes)
chr_state_WGBS = melt(chr_state_WGBS,id.vars=c("chromosome","Sample"))
colnames(chr_state_WGBS)[3:4] = c("State","CpGs")
chr_state_WGBS_total = ddply(chr_state_WGBS,.(chromosome,State),summarise,CpGs=sum(CpGs))

chr_state_DNase = read.table("DNase/chr_DNase.txt",sep=" ")
colnames(chr_state_DNase) = c("Sample","chromosome","Peaks")
chr_state_DNase$chromosome = factor(chr_state_DNase$chromosome,levels=standard_chromosomes)
chr_state_DNase = merge(chr_state_DNase,hg19_genome,by.x="chromosome",by.y="chrom")
chr_state_DNase_total = ddply(chr_state_DNase,.(chromosome),summarise,Peaks=sum(Peaks),Size=sum(as.numeric(size)))

chr_state_H3K27ac = read.table("H3K27ac/chr_H3K27ac.txt",sep=" ")
colnames(chr_state_H3K27ac) = c("Sample","chromosome","Peaks")
chr_state_H3K27ac$chromosome = factor(chr_state_H3K27ac$chromosome,levels=standard_chromosomes)
chr_state_H3K27ac = merge(chr_state_H3K27ac,hg19_genome,by.x="chromosome",by.y="chrom")
chr_state_H3K27ac_total = ddply(chr_state_H3K27ac,.(chromosome),summarise,Peaks=sum(Peaks),Size=sum(as.numeric(size)))

# TE chromosome distribution
# Length of each chromosome
hg19_standard = hg19_genome[which(hg19_genome$chrom %in% standard_chromosomes),]
hg19_standard = ddply(hg19_standard,.(),transform,Prop=size/sum(as.numeric(size)))[2:4]
hg19_standard$chrom = factor(hg19_standard$chrom,levels=standard_chromosomes)

# % of TEs on each chromosome
hg19_TEs = merge(hg19_standard,ddply(rmsk_TE,.(chromosome),summarize,TEs=length(chromosome)),by.x="chrom",by.y="chromosome",all.x=TRUE)
hg19_TEs = ddply(hg19_TEs,.(),transform,Prop_TEs=TEs/NUM_TE)[2:6]

# % of class on each chromosome
test = ddply(rmsk_TE,.(chromosome,class_update),summarise,Class=length(chromosome))
test = ddply(test,.(class_update),transform,Prop_Class=Class/sum(Class))
hg19_class = merge(hg19_standard,test,by.x="chrom",by.y="chromosome",all=TRUE)

test = ddply(melt(rmsk_TE[,c("chromosome",cohorts)],id.var="chromosome"),.(chromosome,variable),summarise,TEs_feature=sum(!is.na(value)))
hg19_features = merge(hg19_TEs,dcast(test,chromosome~variable),by.x="chrom",by.y="chromosome")
hg19_features = melt(hg19_features,id.vars=c("chrom","size","Prop","TEs","Prop_TEs"))
colnames(hg19_features)[6:7] = c("Cohort","TEs_feature")
hg19_features$Prop_TEs_feature = hg19_features$TEs_feature/hg19_features$TEs
hg19_features = split_coding(hg19_features)
hg19_features = ddply(hg19_features,.(Cohort),transform,Zscore=(Prop_TEs_feature-mean(Prop_TEs_feature))/sd(Prop_TEs_feature))

# Chromosome potential
chromosome_potential = ddply(melt(rmsk_TE_measure[,c("chromosome",states)],id.var="chromosome"),.(chromosome,variable),summarise,Samples=mean(na.omit(value)))
colnames(chromosome_potential)[2] = c("State")
chromosome_potential$chromosome = factor(chromosome_potential$chromosome,levels=standard_chromosomes)
chromosome_potential = ddply(chromosome_potential,.(State),transform,Zscore=(Samples-mean(Samples))/sd(Samples))