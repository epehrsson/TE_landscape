# H3K27ac potential
# See 8/3/2017

# Distribution of TEs overlapping DNase peaks
TE_H3K27ac_potential = sample_distribution(TE_H3K27ac_peaks,106,98)
TE_H3K27ac_potential[1,2] = 4430788-1801527
colnames(TE_H3K27ac_potential)[2] = "H3K27ac"

# Cumulative distribution of TEs overlapping DNase peaks
TE_H3K27ac_potential_cum = as.data.frame(cumsum(TE_H3K27ac_potential[which(TE_H3K27ac_potential$Samples %in% seq(1,98)),]$H3K27ac)/sum(TE_H3K27ac_potential[which(TE_H3K27ac_potential$Samples %in% seq(1,98)),]$H3K27ac))
rownames(TE_H3K27ac_potential_cum) = seq(1:98)
TE_H3K27ac_potential_cum$Samples = rownames(TE_H3K27ac_potential_cum)
colnames(TE_H3K27ac_potential_cum)[1] = "Proportion"

# Statistics
TE_H3K27ac_potential_stats = as.data.frame(t(rbind(sum(TE_H3K27ac_potential$H3K27ac[2:99])/sum(TE_H3K27ac_potential$H3K27ac),(sum(as.numeric(TE_H3K27ac_potential$H3K27ac)*seq(0,98))/sum(TE_H3K27ac_potential$H3K27ac))/98,(sum(as.numeric(TE_H3K27ac_potential$H3K27ac[2:99])*seq(1,98))/sum(TE_H3K27ac_potential$H3K27ac))/98)))
colnames(TE_H3K27ac_potential_stats)[1:3] = c("Proportion_ever","Samples_avg_all","Samples_avg_ever")
TE_H3K27ac_potential_stats$State = "yes"

# Proportion of TEs overlapping Dnase peaks, by sample
TE_H3K27ac_peaks_sample = as.data.frame(apply(TE_H3K27ac_peaks[,8:105],2,function(x) sum(x > 0)))
colnames(TE_H3K27ac_peaks_sample) = "H3K27ac"
TE_H3K27ac_peaks_sample$Proportion = TE_H3K27ac_peaks_sample$H3K27ac/4430788
TE_H3K27ac_peaks_sample[which(rownames(TE_H3K27ac_peaks_sample) %in% c("E116","E117","E123","E124","E126","E127")),]$Proportion = TE_H3K27ac_peaks_sample[which(rownames(TE_H3K27ac_peaks_sample) %in% c("E116","E117","E123","E124","E126","E127")),]$H3K27ac/4399208
