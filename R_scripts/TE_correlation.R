# Correlation of TE features with epigenetic measures
# See 4/19/2016, 4/25/2016, 8/24/2016, 8/25/2016, 9/20/2016, 9/21/2016, 9/27/2016, 9/28/2016, 11/4/2016, 11/5/2016, 11/7/2016, 11/18/2016, 12/16/2016, 1/31/2017, 2/1/2017, 2/3/2017, 2/6/2017, 2/9/2017, 2/10/2017, 2/25/2017, 2/27/2017, 2/28/2017, 3/5/2017, 3/8/2017, 5/14/2017, 5/15/2017, 5/16/2017, 5/17/2017, 6/7/2017, 6/14/2017, 6/15/2017, 7/21/2017, 7/24/2017, 8/1/2017, 8/2/2017

# Combine rmsk_TE with other matrices
load("R_scripts/rmsk_TE.RData")
load("R_scripts/chromHMM_TE_state.RData")
load("R_scripts/TE_meth_average.RData")
load("R_scripts/TE_DNase_peaks.RData")
load("R_scripts/TE_H3K27ac_peaks.RData")
load("R_scripts/rna.RData")

rmsk_TE_measure = rmsk_TE
contrasts(rmsk_TE_measure$class_update) <- contr.sum

# Convert feature overlap "NA" to 0 - or categories
rmsk_TE_measure[is.na(rmsk_TE_measure)] = 0
rmsk_TE_measure[,13:32] = apply(rmsk_TE_measure[,13:32],2,function(x) {ifelse(x == 0, "no", "yes")})
rmsk_TE_measure[,13:32] = apply(rmsk_TE_measure[,13:32],2,function(x) x <- factor(x,levels=c("no","yes")))

# Add number of samples in chromHMM state per TE
rmsk_TE_measure = merge(rmsk_TE_measure,chromHMM_TE_state[,c(1:23,25)],by=c("chromosome","start","stop","subfamily","family","class","strand"))

# Add number of samples in methylation state per TE
rmsk_TE_measure = merge(rmsk_TE_measure,TE_meth_average[,c(1:7,45:53)],by=c("chromosome","start","stop","subfamily","family","class","strand"),all.x=TRUE)
rmsk_TE_measure[which(is.na(rmsk_TE_measure$CpGs)),]$CpGs = 0
rmsk_TE_measure$CpGs_per_length = rmsk_TE_measure$CpGs/rmsk_TE_measure$Length

# Add number of samples overlapping DNase peak per TE
rmsk_TE_measure = merge(rmsk_TE_measure,TE_DNase_peaks[,c(1:7,62)],by=c("chromosome","start","stop","subfamily","family","class","strand"),all.x=TRUE)
colnames(rmsk_TE_measure)[length(colnames(rmsk_TE_measure))] = "DNase"
rmsk_TE_measure[which(is.na(rmsk_TE_measure$DNase)),]$DNase = 0

# Add number of samples overlapping H3K27ac peak per TE
rmsk_TE_measure = merge(rmsk_TE_measure,TE_H3K27ac_peaks[,c(1:7,107)],by=c("chromosome","start","stop","subfamily","family","class","strand"),all.x=TRUE)
colnames(rmsk_TE_measure)[length(colnames(rmsk_TE_measure))] = "H3K27ac"
rmsk_TE_measure[which(is.na(rmsk_TE_measure$H3K27ac)),]$H3K27ac = 0

# Add number of samples expressed per TE
rmsk_TE_measure = merge(rmsk_TE_measure,RNA_TE_agnostic[,c(1:7,61:62)],by=c("chromosome","start","stop","subfamily","family","class","strand"))

rmsk_TE_measure = rmsk_TE_measure[,c(1:12,50,59,13:32,33:49,51:58,60:63)]