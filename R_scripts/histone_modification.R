# Loads a dataframe ("histones") with the fold enrichment over input for histone modification ChIP-seq and DHS, 
# plus average DNA methylation level, averaged over all instances of TEs in each chromHMM state, 
# in 50bp bins, extended 5kb from the center of the TEs in each direction
# Excludes the 5_TxWk, 14_ReprPCWk, and 15_Quies states

# Load dataframe of the average level of each epigenetic mark over each 50bp bin, 
# By chromHMM state (all instances of a TE in that state)
histones = read.table("compare_marks/profile_histone/rmsk_TEother_average.txt",sep='\t',
                      col.names=c("Bin","Level","Bases","Mark","State"))
histones$Mark = factor(histones$Mark,levels=levels(histones$Mark)[c(5:6,8,3:4,7,2,1,9)])
histones$State = gsub("8_ZNF.Rpts","8_ZNF/Rpts",histones$State)
histones$State = factor(histones$State,levels=chromHMM_states[c(1:4,6:13)])

# Adjusting scale for DNA methylation level to allow the second y-axis to be visualized opposite the first y-axis
histones[which(histones$Mark == "meth"),]$Level = histones[which(histones$Mark == "meth"),]$Level*12.5

# Converting bin number to bp from the center of the TE for the x-axis
histones$Position = as.numeric(lapply(histones$Bin,function(x) unlist(strsplit(as.character(x),split="Bin_",fixed=TRUE))[2]))
histones$Position = mapvalues(histones$Position,seq(1,200,1),seq(-4975,4975,50))