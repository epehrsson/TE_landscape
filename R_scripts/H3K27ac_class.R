# Class analyses for H3K27ac
# See 7/4/2017

# Contribution
# Contribution of H3K27ac peaks overlapping TEs by class
TE_H3K27ac_class = read.table("H3K27ac/class_H3K27ac_sample.txt",sep='\t')
colnames(TE_H3K27ac_class) = c("class","Sample","Overlap")
TE_H3K27ac_class = dcast(TE_H3K27ac_class,Sample~class)
TE_H3K27ac_class = merge(TE_H3K27ac_class,H3K27ac_stats[,c(1,5)],by=c("Sample"))
TE_H3K27ac_class[,2:9] = t(apply(TE_H3K27ac_class,1,function(x) as.numeric(x[2:9])/as.numeric(x[9])))

# Potential
# Distribution of TEs overlapping H3K27ac peaks, by class
potential_TEother_H3K27ac_class = as.matrix(table(TE_H3K27ac_peaks[,106:107]))

# Adding TEs never overlapping H3K27ac peak
potential_TEother_H3K27ac_class = rbind(rmsk_TEother_class_pi[match(colnames(potential_TEother_H3K27ac_class),rmsk_TEother_class_pi$Class),]$Elements - 
                                          colSums(potential_TEother_H3K27ac_class),potential_TEother_H3K27ac_class)
rownames(potential_TEother_H3K27ac_class) = seq(0,98,1)
potential_TEother_H3K27ac_class = as.data.frame(potential_TEother_H3K27ac_class)
potential_TEother_H3K27ac_class$Total = rowSums(potential_TEother_H3K27ac_class)

# Statistics
potential_TEother_H3K27ac_class_stats = rbind(apply(potential_TEother_H3K27ac_class,2,function(x) sum(x[2:99])/(sum(x)/100)),apply(potential_TEother_H3K27ac_class,2,function(x) sum(x*as.numeric(rownames(potential_TEother_H3K27ac_class)))/sum(x))/0.98,apply(potential_TEother_H3K27ac_class[2:99,],2,function(x) sum(x*as.numeric(rownames(potential_TEother_H3K27ac_class)[2:99]))/sum(x))/0.98)

# Proportion of TEs overlapping DNase peak by class by sample
TE_H3K27ac_peaks_class = aggregate(data=TE_H3K27ac_peaks[,c(8:105,107)],.~class_update,function(x) sum(x > 0))
rownames(TE_H3K27ac_peaks_class) = TE_H3K27ac_peaks_class$class_update
TE_H3K27ac_peaks_class = TE_H3K27ac_peaks_class[,2:99]
TE_H3K27ac_peaks_class = TE_H3K27ac_peaks_class/rmsk_TEother_class_pi[match(rownames(TE_H3K27ac_peaks_class),rmsk_TEother_class_pi$Class),]$Elements