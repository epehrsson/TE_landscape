orthogonal_sets = data.frame(cbind(File=gsub("_vista.txt","",list.files(path = "orthogonal/vista")),
                              State=c(states[c(1:3,6:7)],NA,rep(states[c(1:3,6:7)],11)),
                              Set=c(rep("Constitutive",5),"Inactive",rep("Restricted",5),rep("Shuffled",50)),
                              Iteration=c(rep(NA,11),rep(c(1,10,seq(2,9,1)),each=5))))

vista_val = lapply(list.files(path = "orthogonal/vista",full.names = TRUE),
                   function(x) read.table(x,sep="\t",
                                          col.names=c(TE_coordinates[c(1:4,6,7,5)],"chr_vista","start_vista","stop_vista","element","Validation","Overlap")))
names(vista_val) = as.vector(orthogonal_sets$File)
vista_val = vista_val[sapply(vista_val, function(x) dim(x)[1]) > 0]

vista_val = ldply(vista_val,.id = "File")
vista_val = merge(vista_val,orthogonal_sets,by="File")
  
# Unique VISTA enhancers overlapping TEs in that set
vista_val_fraction = ddply(vista_val,.(Set,State,Iteration),summarise,Total=length(unique(element)),
      Pos=length(unique(element[which(Validation == "positive")])),
      Fraction=Pos/Total)

ggplot(vista_val_fraction,aes(x=Set,y=Fraction,fill=State)) + geom_boxplot() + facet_wrap(~State) + scale_fill_manual(values=all_state_colors,guide=FALSE) + ylim(0,1) + 
  ylab("Proportion of enhancers positively validated") + theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_x_discrete(limits=c("Constitutive","Shuffled","Restricted","Inactive"),labels=c("Constitutively active","Shuffled constitutively active","Restricted activity","Inactive"))
