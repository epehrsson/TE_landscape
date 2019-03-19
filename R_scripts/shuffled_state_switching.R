# Shuffled epigenetic state dynamics

load("R_datasets/shuffled_chromHMM.RData")
load("R_datasets/shuffled_WGBS.RData")

# Total states per TE
shuffled_chromHMM_states = ldply(shuffled_chromHMM_potential,function(x) x[,c("class","States")])
shuffled_chromHMM_states$Iteration = factor(rep(seq(1,10,1),each=NUM_TE))

shuffled_WGBS_average = lapply(seq(1,10,1),function(y) transform(shuffled_WGBS_average[[y]],Iteration = rep(y,unlist(dim(shuffled_WGBS_average[[y]])[1]))))
shuffled_WGBS_states = ldply(shuffled_WGBS_average,function(x) x[,c("class","States","Iteration")])
shuffled_WGBS_states$Iteration = as.factor(shuffled_WGBS_states$Iteration)

## Median total states per TE
table(shuffled_chromHMM_states$States,shuffled_chromHMM_states$Iteration)
ddply(shuffled_chromHMM_states,.(Iteration),summarise,Median=median(States),IQR=IQR(States))
ddply(shuffled_chromHMM_states,.(Iteration,class),summarise,Median=median(States),IQR=IQR(States))

table(shuffled_WGBS_states$States,shuffled_WGBS_states$Iteration)
ddply(shuffled_WGBS_states,.(Iteration),summarise,Median=median(States),IQR=IQR(States))
ddply(shuffled_WGBS_states,.(Iteration,class),summarise,Median=median(States),IQR=IQR(States))

# State switching inter 
## TEs ever in state, chromHMM
shuffled_chromHMM_ever = ldply(shuffled_chromHMM_potential,function(x) apply(x[,8:22],2,function(y) sum(y > 0)))
colnames(shuffled_chromHMM_ever) = chromHMM_states
shuffled_chromHMM_ever$Iteration = seq(1,10,1)
shuffled_chromHMM_ever = melt(shuffled_chromHMM_ever,id.vars="Iteration")
colnames(shuffled_chromHMM_ever)[2:3] = c("State","Count")

## TEs ever in poised promoter states
ldply(shuffled_chromHMM_potential,function(x) dim(x[which(x$X10_TssBiv > 0 | x$X11_BivFlnk > 0),])[1]/NUM_TE)

## Inter state switching
shuffled_chromHMM_inter = lapply(list.files(path="chromHMM/shuffled_TEs",pattern="_inter.txt",full.names=TRUE),function(x) read.table(x,sep='\t',header=TRUE,row.names=1))[c(1,3:10,2)]
shuffled_chromHMM_inter = ldply(shuffled_chromHMM_inter,function(x) melt(as.matrix(x)))
colnames(shuffled_chromHMM_inter) = c("State1","State2","Count")
shuffled_chromHMM_inter$Iteration = rep(seq(1,10,1),each=225)

## Normalize row by number of TEs ever in state
shuffled_chromHMM_inter = merge(shuffled_chromHMM_inter,shuffled_chromHMM_ever,by.x=c("Iteration","State1"),by.y=c("Iteration","State"))
shuffled_chromHMM_inter$Mean_samples = shuffled_chromHMM_inter$Count.x/shuffled_chromHMM_inter$Count.y
shuffled_chromHMM_inter
ddply(shuffled_chromHMM_inter,.(State1,State2),summarise,Proportion=mean(Mean_samples/sample_counts["All","chromHMM"]),SD=sd(Mean_samples/sample_counts["All","chromHMM"]))

ggplot(shuffled_chromHMM_inter,aes(x=State2,y=State1,fill=Mean_samples/sample_counts["All","chromHMM"])) + geom_tile() +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),legend.title = element_blank()) + scale_y_discrete(limits=rev(chromHMM_states)) + 
  scale_fill_gradient(low="white",high="darkmagenta",limits=c(0,1),guide=FALSE) + coord_equal() + facet_wrap(~Iteration)

## Inter state switching, WGBS
shuffled_meth_inter = lapply(list.files(path="WGBS/shuffled",pattern="_inter.txt",full.names=TRUE),function(x) read.table(x,sep='\t',header=TRUE,row.names=1))[c(1,3:10,2)]
shuffled_meth_inter = lapply(shuffled_meth_inter,function(x) transform(x,State1=rownames(x)))
shuffled_meth_inter = ldply(shuffled_meth_inter,function(x) melt(x,id.vars=c("State1","Total")))
colnames(shuffled_meth_inter)[3:4] = c("State2","Count")
shuffled_meth_inter$Iteration = rep(seq(1,10,1),each=16)

## Normalize row by number of TEs ever in state
shuffled_meth_inter$Mean_samples = shuffled_meth_inter$Count/shuffled_meth_inter$Total
shuffled_meth_inter

ggplot(shuffled_meth_inter,aes(x=State2,y=State1,fill=Mean_samples/sample_counts["All","WGBS"])) + geom_tile() +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),legend.title = element_blank()) + scale_y_discrete(limits=rev(meth_states)) + 
  scale_fill_gradient(low="white",high="darkmagenta",limits=c(0,1),guide=FALSE) + coord_equal() + facet_wrap(~Iteration)

# Total states for TEs ever in state
shuffled_chromHMM_dynamics = lapply(shuffled_chromHMM_potential,function(y) apply(y[,8:22],2,function(x) median(y[which(x > 0),]$States)))
shuffled_chromHMM_dynamics = ldply(shuffled_chromHMM_dynamics)
shuffled_chromHMM_dynamics

shuffled_meth_dynamics = lapply(shuffled_WGBS_average,function(y) apply(y[,meth_states],2,function(x) median(y[which(x > 0),]$States)))
shuffled_meth_dynamics = ldply(shuffled_meth_dynamics)
shuffled_meth_dynamics

### Remove shuffled matrices
rm(list=c("shuffled_chromHMM_potential","shuffled_WGBS_average"))
