# Calculates the proportion of TEs in each TE class annotated with each epigenetic state that overlap genic features,
# scaled by the frequency with which the TE is annotated with the state.
# I.e., what proportion of TEs are likely epigenetic passengers?

## state_genic - proportion of TEs in each TE class, annotated with each epigenetic state at each frequency, 
## that overlap RefSeq genes or promoters
## rmsk_TE_genic - proportion of TEs in each class that overlap RefSeq genes or promoters

# Reformat dataframe to list each TE and its length of overlap with introns, exons, and promoters (all and protein-coding only)
# versus the proportion of samples the TE is annotated with each state
rmsk_TE_measure_long = melt(rmsk_TE_measure[,c("class_update","introns","exons","promoters","introns_pc","exons_pc","promoters_pc",states)],
                            id.vars=c("class_update","introns","exons","promoters","introns_pc","exons_pc","promoters_pc"),
                            variable.name="State",value.name="Sample.Proportion")

# By TE class and epigenetic state, count the number of TEs in the state at each frequency (proportion of samples),
# overall, overlapping any RefSeq gene (exon/intron/promoter), or overlapping a protein-coding RefSeq gene
state_genic = ddply(rmsk_TE_measure_long,.(class_update,State,Sample.Proportion),summarise,
            Transcribed=length(class_update),
            All=length(class_update[which(introns > 0 | exons > 0 | promoters > 0)]),
            PC=length(class_update[which(introns_pc > 0 | exons_pc > 0 | promoters_pc > 0)]))
state_genic$Sample.Proportion = as.numeric(state_genic$Sample.Proportion)
## Round the proportion of samples in the state and bin results
state_genic$Bin = round(state_genic$Sample.Proportion,2)
state_genic = ddply(state_genic,.(class_update,State,Bin),summarise,Transcribed=sum(Transcribed),All=sum(All),PC=sum(PC))
## Proportion of TEs overlapping genic features (all or protein-coding only)
state_genic = melt(state_genic,id.vars=c("class_update","State","Bin","Transcribed"),variable.name="Genes",value.name="Genic")
state_genic$Ratio = state_genic$Genic/state_genic$Transcribed

# Proportion of TEs in each class overlapping RefSeq genic features (any or protein-coding only)
rmsk_TE_genic = ddply(rmsk_TE_measure,.(class_update),summarise,
      All=length(class_update[which(introns > 0 | exons > 0 | promoters > 0)])/length(class_update),
      PC=length(class_update[which(introns_pc > 0 | exons_pc > 0 | promoters_pc > 0)])/length(class_update))
rmsk_TE_genic = melt(rmsk_TE_genic,id.var="class_update",variable.name="Genes",value.name="Proportion")

rm(rmsk_TE_measure_long)