# Load shuffled TE potential matrices

# chromHMM
# Number of samples each shuffled TE is in each chromHMM state
shuffled_chromHMM_potential = lapply(list.files(path="chromHMM/shuffled_TEs/",pattern="_potential.txt",full.names = TRUE),function(x) read.table(x,sep='\t',header=TRUE))
shuffled_chromHMM_potential = lapply(shuffled_chromHMM_potential,function(x) transform(x,States = apply(x[,8:22],1,function(x) sum(x > 0))))

# Maximum number of states in a single sample per TE
shuffled_max_intra = lapply(list.files(path="chromHMM/shuffled_TEs/",pattern="_chromHMM_max.txt",full.names = TRUE),function(x) read.table(x,sep='\t'))
shuffled_max_intra = lapply(shuffled_max_intra, setNames, nm =c("chromosome","start","stop","subfamily","class","family","strand","Max_states_intra"))
for (i in 1:10){
  shuffled_chromHMM_potential[[i]] = merge(shuffled_chromHMM_potential[[i]],shuffled_max_intra[[i]],by=c("chromosome","start","stop","subfamily","class","family","strand"))
}
rm(shuffled_max_intra)

# WGBS
# Methylation level for shuffled TEs
shuffled_WGBS_average = lapply(list.files(path="WGBS/shuffled/",pattern="_Meth.bed",full.names = TRUE),function(x) read.table(x,sep='\t',header=TRUE))

# Number of CpGs per TE
shuffled_WGBS_CpG = lapply(list.files(path="WGBS/shuffled/",pattern="TE_CpG_count_",full.names = TRUE),read.table(x,sep='\t'))
shuffled_WGBS_CpG = lapply(shuffled_WGBS_CpG, setNames, nm =c("chromosome","start","stop","subfamily","class","family","strand","CpGs"))

# Methylation level for TEs with at least one CpG
for (i in 1:10){
  shuffled_WGBS_average[[i]] = merge(shuffled_WGBS_average[[i]],shuffled_WGBS_CpG[[i]],by=c("chromosome","start","stop","subfamily","class","family","strand"))
}
rm(shuffled_WGBS_CpG)

# Number of samples TE is in each methylation state (IMR90, no IMR90)
shuffled_WGBS_average = lapply(shuffled_WGBS_average,function(x) transform(x,Hypomethylated = apply(x,1,function(x) sum(x[8:44] < 0.3,na.rm=TRUE)),Hypermethylated = apply(x,1,function(x) sum(x[8:44] > 0.7,na.rm=TRUE)),Intermediate = apply(x,1,function(x) sum(x[8:44] <= 0.7 & x[8:44] >= 0.3,na.rm=TRUE)),Missing = apply(x,1,function(x) sum(is.na(x[8:44]))),Hypomethylated_noIMR90 = apply(x,1,function(x) sum(x[c(8:17,19:44)] < 0.3,na.rm=TRUE)),Hypermethylated_noIMR90 = apply(x,1,function(x) sum(x[c(8:17,19:44)] > 0.7,na.rm=TRUE)),Intermediate_noIMR90 = apply(x,1,function(x) sum(x[c(8:17,19:44)] <= 0.7 & x[c(8:17,19:44)] >= 0.3,na.rm=TRUE)),Missing_noIMR90 = apply(x,1,function(x) sum(is.na(x[c(8:17,19:44)])))))

# DNase
# Number of overlaps between shuffled TEs, DNase peaks
shuffled_DNase_potential = lapply(list.files(path="DNase/shuffled/",pattern="_DNase_peaks.txt",full.names = TRUE),function(x) read.table(x,sep='\t'))
shuffled_DNase_potential = lapply(shuffled_DNase_potential, setNames, nm =c("chromosome","start","stop","subfamily","class","family","strand","Sample","Peaks","Overlap"))
shuffled_DNase_potential = lapply(shuffled_DNase_potential,function(x) dcast(shuffled_DNase_potential,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks"))
shuffled_DNase_potential = lapply(shuffled_DNase_potential,function(x) x[is.na(x)] <- 0)                                  

# Number of samples a TE overlaps a DNase peak
shuffled_DNase_potential = lapply(shuffled_DNase_potential,function(x) transform(x,Samples = apply(x,1,function(x) sum(as.numeric(x[8:60]) > 0))))

# H3K27ac
# Number of overlaps between shuffled TEs, H3K27ac peaks
shuffled_H3K27ac_potential = lapply(list.files(path="H3K27ac/shuffled/",pattern="_H3K27ac_peaks.txt",full.names = TRUE),function(x) read.table(x,sep='\t'))
shuffled_H3K27ac_potential = lapply(shuffled_H3K27ac_potential, setNames, nm =c("chromosome","start","stop","subfamily","class","family","strand","Sample","Peaks","Overlap"))
shuffled_H3K27ac_potential = lapply(shuffled_H3K27ac_potential,function(x) dcast(shuffled_H3K27ac_potential,chromosome+start+stop+subfamily+family+class+strand~Sample,value.var="Peaks"))
shuffled_H3K27ac_potential = lapply(shuffled_H3K27ac_potential,function(x) x[is.na(x)] <- 0)                                  

# Number of samples a TE overlaps a H3K27ac peak
shuffled_H3K27ac_potential = lapply(shuffled_H3K27ac_potential,function(x) transform(x,Samples = apply(x,1,function(x) sum(as.numeric(x[8:105]) > 0))))

save(list=c("shuffled_chromHMM_potential","shuffled_WGBS_average","shuffled_DNase_potential","shuffled_H3K27ac_potential"),file="R_scripts/shuffled_matrices.RData")