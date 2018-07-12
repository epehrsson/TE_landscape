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
enrich_H3K27ac = merge(enrich_H3K27ac,metadata[,c("Sample",sample_categories)])

## Select training and validation sets (70% and 30%)
train <- sample(nrow(enrich_H3K27ac), 0.7*nrow(enrich_H3K27ac), replace = FALSE)
TrainSet <- droplevels(enrich_H3K27ac[train,])
ValidSet <- droplevels(enrich_H3K27ac[-train,])

# Number of variables to use
result <- rfcv(enrich_H3K27ac[,c(2:936)],enrich_H3K27ac$Germline,cv.fold=5)

# Create a Random Forest model with default parameters
model_H3K27ac_Group <- randomForest(Group ~ ., data = TrainSet[,c(2:937)], importance = TRUE, mtry = 120)
model_H3K27ac_Group

# Predicting on validation set
predValid <- predict(model_H3K27ac_Group, ValidSet, type = "class")

# Checking classification accuracy
mean(predValid == as.vector(ValidSet$Group))
table(predValid, ValidSet$Group) 

# To check important variables
importance(model_H3K27ac_Group)        
varImpPlot(model_H3K27ac_Group)
