# Clustering by sample
sample_cluster_3D = dcast(subfamily_state_sample_combined[,c("subfamily","State","Sample","Enrichment")],
                          subfamily+State~Sample,value.var = "Enrichment")
rownames(sample_cluster_3D) = paste(sample_cluster_3D$subfamily,sample_cluster_3D$State,sep="_")
sample_cluster_3D = sample_cluster_3D[,3:129]
sample_cluster = hclust(dist(t(2^sample_cluster_3D)))

# Clustering by subfamily
subfamily_cluster_3D = dcast(subfamily_state_sample_combined[,c("subfamily","State","Sample","Enrichment")],
                             Sample+State~subfamily,value.var = "Enrichment")
rownames(subfamily_cluster_3D) = paste(subfamily_cluster_3D$Sample,subfamily_cluster_3D$State,sep="_")
subfamily_cluster_3D = subfamily_cluster_3D[,3:970]
subfamily_cluster = hclust(dist(t(2^subfamily_cluster_3D)))

# Clustering by state
state_cluster_3D = dcast(subfamily_state_sample_combined[,c("subfamily","State","Sample","Enrichment")],
                         Sample+subfamily~State,value.var = "Enrichment")
rownames(state_cluster_3D) = paste(state_cluster_3D$Sample,state_cluster_3D$subfamily,sep="_")
state_cluster_3D = state_cluster_3D[,3:23]
state_cluster = hclust(dist(t(2^state_cluster_3D)))

test = subfamily_state_sample_combined
test$State = factor(test$State,levels=levels(subfamily_state_sample_combined$State)[state_cluster$order])
test$subfamily = factor(test$subfamily,levels(subfamily_state_sample_combined$subfamily)[subfamily_cluster$order])
test$Sample = factor(test$Sample,levels(subfamily_state_sample_combined$Sample)[sample_cluster$order])
plot_ly(test, x = ~subfamily, y = ~State, z = ~Sample,
        marker = list(size = ~Enrichment, showscale = TRUE)) %>% add_markers() %>% 
  layout(scene = list(xaxis = list(type="category"),yaxis = list(type="category",title="Subfamily"),zaxis = list(type="category")))
