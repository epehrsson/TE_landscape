tf = read.table("TF_matrix.txt",sep='\t',header=TRUE)
tf = tf[,1:4]
tf$Target_PC = as.numeric(gsub("%","",tf$Target_PC))
tf$Background_PC = as.numeric(gsub("%","",tf$Background_PC))
tf = melt(tf,id.vars=c("subfamily","TF"))

ggplot(tf,aes(x=TF,y=value,fill=variable)) + geom_bar(stat="identity",position="dodge") + facet_grid(~subfamily,scales="free_x",space="free_x") + ylab("% TEs with motif") + theme(axis.text.x = element_text(angle=90,hjust=1),axis.title.x = element_blank()) + scale_fill_discrete(guide=FALSE)
