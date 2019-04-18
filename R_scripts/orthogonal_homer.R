homer_known = lapply(list.files(path="orthogonal/homer/",pattern="knownResults.txt",recursive=TRUE,full.names=TRUE),function(x)
                     read.table(x,sep='\t',header=TRUE,comment.char=""))
names(homer_known) = list.files(path="orthogonal/homer/",pattern="knownResults.txt",recursive=TRUE,full.names=TRUE)
homer_known = lapply(homer_known,setNames,nm=c("Motif Name","Consensus","P-value","Log P-value","q-value (Benjamini)",
                                               "Target Seqs with Motif","% of Target Sequences with Motif",
                                               "Background Seqs with Motif","% of Background Sequences with Motif"))

homer_known = ldply(homer_known,.id = "File")
homer_known$File = gsub("/knownResults.txt","",gsub("orthogonal/homer//","",homer_known$File))
homer_known = merge(homer_known,orthogonal_sets,by="File")

homer_known$`% of Target Sequences with Motif` = as.numeric(gsub("%","",homer_known$`% of Target Sequences with Motif`))
homer_known$`% of Background Sequences with Motif` = as.numeric(gsub("%","",homer_known$`% of Background Sequences with Motif`))

# Proportion with TATA box motif
table(homer_known[grep("TATA",homer_known$`Motif Name`),]$`q-value (Benjamini)`)
ggplot(homer_known[grep("TATA",homer_known$`Motif Name`),],aes(x=Set,y=`% of Target Sequences with Motif`,fill=State)) + geom_boxplot() + 
  facet_wrap(~State) + ylab("% of seqs with motif") + ylim(0,100) + scale_fill_manual(values=all_state_colors,guide=FALSE) +
  scale_x_discrete(limits=c("Constitutive","Shuffled","Restricted","Inactive"),labels=c("Constitutively active","Shuffled constitutively active","Restricted activity","Inactive")) + 
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))