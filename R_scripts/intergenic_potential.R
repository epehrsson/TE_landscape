rmsk_TE_measure_long = melt(rmsk_TE_measure[,c("class_update","introns","exons","promoters","introns_pc","exons_pc","promoters_pc",states)],
                            id.vars=c("class_update","introns","exons","promoters","introns_pc","exons_pc","promoters_pc"))
colnames(rmsk_TE_measure_long)[8:9] = c("State","Sample.Proportion")

state_genic = ddply(rmsk_TE_measure_long,.(class_update,State,Sample.Proportion),summarise,
            Transcribed=length(class_update),
            All=length(class_update[which(introns > 0 | exons > 0 | promoters > 0)]),
            PC=length(class_update[which(introns_pc > 0 | exons_pc > 0 | promoters_pc > 0)]))
state_genic$Sample.Proportion = as.numeric(state_genic$Sample.Proportion)
state_genic$Bin = round(state_genic$Sample.Proportion,2)
state_genic = ddply(state_genic,.(class_update,State,Bin),summarise,
                           Transcribed=sum(Transcribed),All=sum(All),PC=sum(PC))
state_genic = melt(state_genic,id.vars=c("class_update","State","Bin","Transcribed"))
colnames(state_genic)[5:6] = c("Genes","Genic")
state_genic$Ratio = state_genic$Genic/state_genic$Transcribed

rmsk_TE_genic = ddply(rmsk_TE_measure,.(class_update),summarise,
      All=length(class_update[which(introns > 0 | exons > 0 | promoters > 0)])/length(class_update),
      PC=length(class_update[which(introns_pc > 0 | exons_pc > 0 | promoters_pc > 0)])/length(class_update))
rmsk_TE_genic = melt(rmsk_TE_genic,id.var="class_update")
colnames(rmsk_TE_genic)[2:3] = c("Genes","Proportion")