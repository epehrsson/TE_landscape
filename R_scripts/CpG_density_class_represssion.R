# Quasi-poisson models of the number of samples a TE is annotated with the 9_Het state or hypermethylated
# Versus class (LTR/SINE) and CpG density

## measure_het - Number of samples each TE is annotated with the 9_Het state and CpG density
## model_het1 - Quasi-poisson model of the impact of class on samples annotated with 9_Het
## model_het2 - Quasi-poisson model of the impact of class and CpG density on samples annotated with 9_Het
## measure_hyper - Number of samples each TE is hypermethylated and CpG density
## model_hyper1 - Quasi-poisson model of the impact of class on samples hypermethylated
## model_hyper2 - Quasi-poisson model of the impact of class and CpG density on samples hypermethylated

# 9_Het
## Number of samples each TE is annotated with the 9_Het state and CpG density (CpGs/kbp)
measure_het = merge(rmsk_TE[,c(TE_coordinates,"class_update","CpGs_per_length")],
                    chromHMM_TE_state[,c(TE_coordinates,"class_update","9_Het")],by=c(TE_coordinates,"class_update"))

## Quasi-poisson model of the impact of class (LTR/SINE elements only) on samples annotated with 9_Het
model_het1 = glm(`9_Het`~class_update,data=measure_het[which(measure_het$class_update %in% c("LTR","SINE")),],family = quasipoisson())

## Quasi-poisson model of the impact of class (LTR/SINE elements only) and CpG density on samples annotated with 9_Het
model_het2 = glm(`9_Het`~class_update+I(CpGs_per_length*1000),data=measure_het[which(measure_het$class_update %in% c("LTR","SINE")),],
                 family = quasipoisson())

# Hypermethylation
## Number of samples each TE is hypermethylated and CpG density (CpGs/kbp)
measure_hyper = merge(rmsk_TE[,c(TE_coordinates,"class_update","CpGs_per_length")],
                    TE_meth_average[,c(TE_coordinates,"class_update","Hypermethylated")],by=c(TE_coordinates,"class_update"))

## Quasi-poisson model of the impact of class (LTR/SINE elements only) on samples hypermethylated
model_hyper1 = glm(Hypermethylated~class_update,data=measure_hyper[which(measure_hyper$class_update %in% c("LTR","SINE")),],family = quasipoisson())

## Quasi-poisson model of the impact of class (LTR/SINE elements only) and CpG density on samples hypermethylated
model_hyper2 = glm(Hypermethylated~class_update+I(CpGs_per_length*1000),
                   data=measure_hyper[which(measure_hyper$class_update %in% c("LTR","SINE")),],family = quasipoisson())