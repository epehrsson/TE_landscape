# Correlation of TE features with epigenetic measures
# See 4/19/2016, 4/25/2016, 8/24/2016, 8/25/2016, 9/20/2016, 9/21/2016, 9/27/2016, 9/28/2016, 11/4/2016, 11/5/2016, 11/7/2016, 11/18/2016, 12/16/2016, 1/31/2017, 2/1/2017, 
# 2/3/2017, 2/6/2017, 2/9/2017, 2/10/2017, 2/25/2017, 2/27/2017, 2/28/2017, 3/5/2017, 3/8/2017, 5/14/2017, 5/15/2017, 5/16/2017, 5/17/2017, 6/7/2017, 6/14/2017, 6/15/2017, 
# 7/21/2017, 7/24/2017, 8/1/2017, 8/2/2017

# Combine rmsk_TE with other matrices
#load("R_datasets/rmsk_TE.RData")
#load("R_datasets/chromHMM_TE_state.RData")
#load("R_datasets/TE_meth_average.RData")
#load("R_datasets/TE_DNase_peaks.RData")
#load("R_datasets/TE_H3K27ac_peaks.RData")
#load("R_datasets/rna.RData")

rmsk_TE_measure = rmsk_TE
contrasts(rmsk_TE_measure$class_update) <- contr.sum
rmsk_TE_measure[is.na(rmsk_TE_measure)] = 0

# Add number of samples in chromHMM state per TE
rmsk_TE_measure = merge(rmsk_TE_measure,chromHMM_TE_state[,c(TE_coordinates,chromHMM_states_X,"States","Max_states_intra")],by=TE_coordinates)

# Add number of samples in methylation state per TE
rmsk_TE_measure = merge(rmsk_TE_measure,TE_meth_average[,c(TE_coordinates,"CpGs",meth_states)],by=TE_coordinates,all.x=TRUE)
rmsk_TE_measure[which(is.na(rmsk_TE_measure$CpGs)),]$CpGs = 0
rmsk_TE_measure$CpGs_per_length = rmsk_TE_measure$CpGs/rmsk_TE_measure$Length

# Add number of samples overlapping DNase peak per TE
rmsk_TE_measure = merge(rmsk_TE_measure,TE_DNase_peaks[,c(TE_coordinates,"Samples")],by=TE_coordinates,all.x=TRUE)
colnames(rmsk_TE_measure)[length(colnames(rmsk_TE_measure))] = "DNase"
rmsk_TE_measure[which(is.na(rmsk_TE_measure$DNase)),]$DNase = 0

# Add number of samples overlapping H3K27ac peak per TE
rmsk_TE_measure = merge(rmsk_TE_measure,TE_H3K27ac_peaks[,c(TE_coordinates,"Samples")],by=TE_coordinates,all.x=TRUE)
colnames(rmsk_TE_measure)[length(colnames(rmsk_TE_measure))] = "H3K27ac"
rmsk_TE_measure[which(is.na(rmsk_TE_measure$H3K27ac)),]$H3K27ac = 0

# Add number of samples expressed per TE
rmsk_TE_measure = merge(rmsk_TE_measure,RNA_TE_agnostic[,c(TE_coordinates,"Expressed_samples","Max_expression")],by=TE_coordinates)

rmsk_TE_measure = rmsk_TE_measure[,c(TE_coordinates,"class_update",measure_metrics,cohorts,measure_states,measure_states_extra)]
