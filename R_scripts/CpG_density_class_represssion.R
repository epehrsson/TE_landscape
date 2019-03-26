# 9_Het
measure_het = merge(rmsk_TE[,c(TE_coordinates,"class_update","CpGs_per_length")],
                    chromHMM_TE_state[,c(TE_coordinates,"class_update","9_Het")],by=c(TE_coordinates,"class_update"))

model_het1 = glm(`9_Het`~class_update,data=measure_het[which(measure_het$class_update %in% c("LTR","SINE")),],family = quasipoisson())
model_het2 = glm(`9_Het`~class_update+I(CpGs_per_length*1000),data=measure_het[which(measure_het$class_update %in% c("LTR","SINE")),],
                 family = quasipoisson())
model_het3 = glm(`9_Het`~class_update+I(CpGs_per_length*1000)+class_update*I(CpGs_per_length*1000),
                data=measure_het[which(measure_het$class_update %in% c("LTR","SINE")),],family = quasipoisson())

# Hyper
measure_hyper = merge(rmsk_TE[,c(TE_coordinates,"class_update","CpGs_per_length")],
                    TE_meth_average[,c(TE_coordinates,"class_update","Hypermethylated")],by=c(TE_coordinates,"class_update"))

model_hyper1 = glm(Hypermethylated~class_update,data=measure_hyper[which(measure_hyper$class_update %in% c("LTR","SINE")),],family = quasipoisson())
model_hyper2 = glm(Hypermethylated~class_update+I(CpGs_per_length*1000),
                   data=measure_hyper[which(measure_hyper$class_update %in% c("LTR","SINE")),],family = quasipoisson())
model_hyper3 = glm(Hypermethylated~class_update+I(CpGs_per_length*1000)+class_update*I(CpGs_per_length*1000),
                data=measure_hyper[which(measure_hyper$class_update %in% c("LTR","SINE")),],family = quasipoisson())
