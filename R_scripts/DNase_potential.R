# DNase potential
# See 8/3/2017

# Distribution of TEs overlapping DNase peaks (needs matrix)
TE_DNase_potential = sample_distribution(TE_DNase_peaks,61,53)
TE_DNase_potential[1,2] = 4430788-1840984
colnames(TE_DNase_potential)[2] = "DNase"

# Cumulative distribution of TEs overlapping DNase peaks
TE_DNase_potential_cum = as.data.frame(cumsum(TE_DNase_potential[which(TE_DNase_potential$Samples %in% seq(1,53)),]$DNase)/sum(TE_DNase_potential[which(TE_DNase_potential$Samples %in% seq(1,53)),]$DNase))
rownames(TE_DNase_potential_cum) = seq(1:53)
TE_DNase_potential_cum$Samples = rownames(TE_DNase_potential_cum)
colnames(TE_DNase_potential_cum)[1] = "Proportion"

# Statistics
TE_DNase_potential_stats = as.data.frame(t(rbind(sum(TE_DNase_potential$DNase[2:54])/sum(TE_DNase_potential$DNase),(sum(as.numeric(TE_DNase_potential$DNase)*seq(0,53))/sum(TE_DNase_potential$DNase))/53,(sum(as.numeric(TE_DNase_potential$DNase[2:54])*seq(1,53))/sum(TE_DNase_potential$DNase))/53)))
colnames(TE_DNase_potential_stats)[1:3] = c("Proportion_ever","Samples_avg_all","Samples_avg_ever")
TE_DNase_potential_stats$State = "yes"

# Proportion of TEs overlapping Dnase peaks, by sample
TE_DNase_peaks_sample = as.data.frame(apply(TE_DNase_peaks[,8:60],2,function(x) sum(x > 0)))
colnames(TE_DNase_peaks_sample) = "DNase"
TE_DNase_peaks_sample$Proportion = TE_DNase_peaks_sample$DNase/4430788
TE_DNase_peaks_sample[which(rownames(TE_DNase_peaks_sample) %in% c("E116","E117","E123","E124","E126","E127")),]$Proportion = TE_DNase_peaks_sample[which(rownames(TE_DNase_peaks_sample) %in% c("E116","E117","E123","E124","E126","E127")),]$DNase/4399208
