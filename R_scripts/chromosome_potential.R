# Chromosome potential

# Overall epigenetic state by chromosome
chrom_state = read.table("chromHMM/chromosome_states.txt",sep='\t')
colnames(chrom_state) = c("Sample","chromosome","State","Bases")
chrom_state$chromosome = factor(chrom_state$chromosome,levels=standard_chromosomes)
chrom_state$State = factor(chrom_state$State,levels=chromHMM_states)
chrom_state_total = ddply(chrom_state,.(chromosome,State),summarise,Bases=sum(as.numeric(Bases)))

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