# Random Forests
## Matrix of enrichments
enrich_H3K27ac = dcast(subfamily_state_sample_combined[which(subfamily_state_sample_combined$State == "H3K27ac" &
                                                               subfamily_state_sample_combined$Count > THRESHOLD_IK_MEMBER),
                                                       c("subfamily","Sample","Enrichment")],Sample~subfamily,value.var = "Enrichment")

## Remove subfamilies under member threshold in some samples
enrich_H3K27ac = enrich_H3K27ac[,which(count_na(enrich_H3K27ac) == 0),]

## Fix names 
colnames(enrich_H3K27ac) = make.names(colnames(enrich_H3K27ac))

## To remove -Inf
enrich_H3K27ac[,2:ncol(enrich_H3K27ac)] = 2^enrich_H3K27ac[,2:ncol(enrich_H3K27ac)]

## Include metadata classifications
enrich_H3K27ac = droplevels(merge(enrich_H3K27ac,metadata[,c("Sample",sample_categories)]))

# Number of variables to use
result <- rfcv(trainx=enrich_H3K27ac[,c(2:936)],trainy=enrich_H3K27ac$Type,cv.fold=5)
result$error.cv

# Create a Random Forest model with default parameters
model_H3K27ac_Group <- randomForest(Group ~ ., data = enrich_H3K27ac[,c(2:937)], importance = TRUE, mtry = 150, ntree=500)
model_H3K27ac_Group

# To check important variables
importance(model_H3K27ac_Group)        
varImpPlot(model_H3K27ac_Group)
