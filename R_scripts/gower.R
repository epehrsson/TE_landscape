c[which(c$p.value < 0.005),]

ggplot(test[which(test$State == "13_ReprPC" & test$Category == "Type"),],
       aes(x=reorder(Grouping,Proportion,median),y=Proportion,fill=Grouping)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle=90)) + scale_fill_manual(values=type_colors)

# Load library
library(cluster)

# A very small subfamily (28 members) with only one state per sample per element
ucon21 = read.table("all_UCON21.txt")
colnames(ucon21) = c(TE_coordinates[c(1:4,6,5,7)],"Sample","Overlap","State","Category")
ucon21 = dcast(ucon21,chromosome+start+stop+subfamily+class+family+strand~Sample,value.var="State")
ucon21[,8:134] = lapply(ucon21[,8:134], function(x) factor(x,levels=chromHMM_states))

# Get a distance matrix using the Gower distance metric for categorical data
# https://medium.com/@anastasia.reusova/hierarchical-clustering-on-categorical-data-in-r-a27e578f2995
# https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/daisy.html
a = daisy(ucon21[,8:134],metric="gower")
y = daisy(as.data.frame(t(ucon21[,8:134])),metric="gower")

# Get clustering from the distance matrix
b = hclust(a, method="complete")
z = hclust(y, method="complete")

# Plotting matrix by replacing states with numbers - map back to colors
heatmap(t(ldply(ucon21[,8:134],function(x) mapvalues(x,chromHMM_states,seq(1,15,1)))[,2:29]),Rowv=as.dendrogram(b),Colv=as.dendrogram(z),
        labCol=metadata$Sample,ColSideColors=group_colors[as.vector(metadata$Group)],col=as.vector(chromHMM_colors),breaks=seq(1,16,1))
