# Class analyses with DNase
# See 6/5/2017

# Potential
# Distribution of TEs overlapping DNase peaks, by class
potential_TEother_DNase_class = as.matrix(table(TE_DNase_peaks[,61:62]))

# Adding TEs never overlapping DNase peak
potential_TEother_DNase_class = rbind(rmsk_TEother_class_pi[match(colnames(potential_TEother_DNase_class),rmsk_TEother_class_pi$Class),]$Elements - colSums(potential_TEother_DNase_class),potential_TEother_DNase_class)
rownames(potential_TEother_DNase_class) = seq(0,53,1)
potential_TEother_DNase_class = as.data.frame(potential_TEother_DNase_class)
potential_TEother_DNase_class$Total = rowSums(potential_TEother_DNase_class)

# Statistics
potential_TEother_DNase_class_stats = rbind(apply(potential_TEother_DNase_class,2,function(x) sum(x[2:54])/(sum(x)/100)),apply(potential_TEother_DNase_class,2,function(x) sum(x*as.numeric(rownames(potential_TEother_DNase_class)))/sum(x))/0.53,apply(potential_TEother_DNase_class[2:54,],2,function(x) sum(x*as.numeric(rownames(potential_TEother_DNase_class)[2:54]))/sum(x))/0.53)

# Contribution
# Contribution of DNase peaks overlapping TEs by class
TE_DNase_class = read.table("DNase_peaks/class_DNase_sample.txt",sep='\t')
colnames(TE_DNase_class) = c("class","Sample","Overlap")
TE_DNase_class = dcast(TE_DNase_class,Sample~class)
TE_DNase_class = merge(TE_DNase_class,DNase_stats[,c(1,5)],by=c("Sample"))
TE_DNase_class[,2:9] = t(apply(TE_DNase_class,1,function(x) as.numeric(x[2:9])/as.numeric(x[9])))

# Proportion of TEs overlapping DNase peak by class by sample
TE_DNase_peaks_class = aggregate(data=TE_DNase_peaks[,c(8:60,62)],.~class_update,function(x) sum(x > 0))
rownames(TE_DNase_peaks_class) = TE_DNase_peaks_class$class_update
TE_DNase_peaks_class = TE_DNase_peaks_class[,2:54]
TE_DNase_peaks_class = TE_DNase_peaks_class/rmsk_TEother_class_pi[match(rownames(TE_DNase_peaks_class),rmsk_TEother_class_pi$Class),]$Elements