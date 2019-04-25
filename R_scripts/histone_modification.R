# Histone modification ChIP-seq, DNase, and methylation over TEs in each chromHMM state
# See 12/13/2016, 3/8/2017, 8/7/2017

histones = read.table("compare_marks/profile_histone/rmsk_TEother_average.txt",sep='\t')
colnames(histones) = c("Bin","Level","Bases","Mark","State")
histones$Mark = factor(histones$Mark,levels=levels(histones$Mark)[c(5:6,8,3:4,7,2,1,9)])
histones$State = gsub("8_ZNF.Rpts","8_ZNF/Rpts",histones$State)
histones$State = factor(histones$State,levels=chromHMM_states[c(1:4,6:13)])

# Adjusting DNA methylation level for the second axis
histones[which(histones$Mark == "meth"),]$Level = histones[which(histones$Mark == "meth"),]$Level*12.5

# Adding x-axis
histones$Position = as.numeric(lapply(histones$Bin,function(x) unlist(strsplit(as.character(x),split="Bin_",fixed=TRUE))[2]))
histones$Position = mapvalues(histones$Position,seq(1,200,1),seq(-4975,4975,50))