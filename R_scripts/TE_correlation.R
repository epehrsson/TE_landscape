# Creates a single data frame ("rmsk_TE_measure") with individual TE characteristics and
# the proportion of samples each TE is annotated with each epigenetic state

# Dataframe of individual TEs with characteristics such as age, mappability, and genic feature overlap
rmsk_TE_measure = rmsk_TE

## Make class levels equivalent
contrasts(rmsk_TE_measure$class_update) <- contr.sum

## Convert feature overlap from NA to 0
rmsk_TE_measure[is.na(rmsk_TE_measure)] = 0

# Add the proportion of samples each TE is annotated with each chromHMM state
# Plus the total number of unique chromHMM states each TE is annotated with across all samples
# And the maximum number of chromHMM states each TE is annotated with in a single sample
rmsk_TE_measure = merge(rmsk_TE_measure,chromHMM_TE_state[,c(TE_coordinates,chromHMM_states,"States","Max_states_intra")],by=TE_coordinates)
rmsk_TE_measure = rename(rmsk_TE_measure,c("States"="States.chromHMM"))
rmsk_TE_measure[,chromHMM_states] = apply(rmsk_TE_measure[,chromHMM_states],2,function(x) as.numeric(ifelse(rmsk_TE_measure$chromosome != "chrY",x/sample_counts["All","chromHMM"],x/sample_counts["chrY","chromHMM"])))

# Add the proportion of samples each TE is annotated with each methylation state
# Plus the total number of unique methylation states each TE is annotated with across all samples
rmsk_TE_measure = merge(rmsk_TE_measure,TE_meth_average[,c(TE_coordinates,meth_states,"States")],by=TE_coordinates,all.x=TRUE)
rmsk_TE_measure = rename(rmsk_TE_measure,c("States"="States.WGBS"))
rmsk_TE_measure[,meth_states] = apply(rmsk_TE_measure[,meth_states],2,function(x) as.numeric(x/sample_counts["All","WGBS"]))

# Add the proportion of samples each TE overlaps the summit of a DHS peak
rmsk_TE_measure = merge(rmsk_TE_measure,TE_DNase_peaks[,c(TE_coordinates,"Samples")],by=TE_coordinates,all.x=TRUE)
colnames(rmsk_TE_measure)[length(colnames(rmsk_TE_measure))] = "DNase"
rmsk_TE_measure[which(is.na(rmsk_TE_measure$DNase)),]$DNase = 0
rmsk_TE_measure$DNase = as.numeric(ifelse(rmsk_TE_measure$chromosome != "chrY",rmsk_TE_measure$DNase/sample_counts["All","DNase"],rmsk_TE_measure$DNase/sample_counts["chrY","DNase"]))

# Add the proportion of samples each TE overlaps the summit of an H3K27ac peak
rmsk_TE_measure = merge(rmsk_TE_measure,TE_H3K27ac_peaks[,c(TE_coordinates,"Samples")],by=TE_coordinates,all.x=TRUE)
colnames(rmsk_TE_measure)[length(colnames(rmsk_TE_measure))] = "H3K27ac"
rmsk_TE_measure[which(is.na(rmsk_TE_measure$H3K27ac)),]$H3K27ac = 0
rmsk_TE_measure$H3K27ac = as.numeric(ifelse(rmsk_TE_measure$chromosome != "chrY",rmsk_TE_measure$H3K27ac/sample_counts["All","H3K27ac"],rmsk_TE_measure$H3K27ac/sample_counts["chrY","H3K27ac"]))

# Add the proportion of samples each TE is expressed RPKM > 1
# Plus the maximum expression across all samples
rmsk_TE_measure = merge(rmsk_TE_measure,RNA_TE[,c(TE_coordinates,"Expressed_samples","Max_expression")],by=TE_coordinates)
rmsk_TE_measure$Expressed_samples = as.numeric(ifelse(rmsk_TE_measure$chromosome != "chrY",rmsk_TE_measure$Expressed_samples/sample_counts["All","RNA"],rmsk_TE_measure$Expressed_samples/sample_counts["chrY","RNA"]))

# Rearrange the dataframe and restrict TE characteristics to length, mappability, JC distance from subfamily consensus, number of CpGs, and CpG density
rmsk_TE_measure = rmsk_TE_measure[,c(TE_coordinates,"class_update",measure_metrics,cohorts,states,measure_states_extra)]