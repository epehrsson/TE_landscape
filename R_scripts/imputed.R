imputed = read.table("/scratch/ecp/TE_landscape/imputed/counts_imputed.txt",sep="\t",quote="")
colnames(imputed) = c("State","Imputed.State","Count")
imputed$State = factor(imputed$State,levels=chromHMM_states)
imputed$Imputed.State = factor(imputed$Imputed.State,levels=levels(imputed$Imputed.State)[c(1,12,19:25,2:11,13:18)])

imputed_states = setNames(c(rgb(255,0,0,maxColorValue=255),rgb(255,69,0,maxColorValue=255),rgb(255,69,0,maxColorValue=255),
                            rgb(255,69,0,maxColorValue=255),rgb(0,128,0,maxColorValue=255),rgb(0,128,0,maxColorValue=255),rgb(0,128,0,maxColorValue=255),
                            rgb(0,150,0,maxColorValue=255),rgb(194,225,5,maxColorValue=255),rgb(194,225,5,maxColorValue=255),rgb(194,225,5,maxColorValue=255),
                            rgb(194,225,5,maxColorValue=255),rgb(255,195,77,maxColorValue=255),rgb(255,195,77,maxColorValue=255),rgb(255,195,77,maxColorValue=255),
                            rgb(255,255,0,maxColorValue=255),rgb(255,255,0,maxColorValue=255),rgb(255,255,0,maxColorValue=255),rgb(255,255,102,maxColorValue=255),
                            rgb(102,205,170,maxColorValue=255),rgb(138,145,208,maxColorValue=255),rgb(230,184,183,maxColorValue=255),
                            rgb(112,48,160,maxColorValue=255),rgb(128,128,128,maxColorValue=255),rgb(255,255,255,maxColorValue=255)),
                          c("1_TssA","2_PromU","3_PromD1","4_PromD2","5_Tx5'","6_Tx","7_Tx3'","8_TxWk",
                            "9_TxReg","10_TxEnh5'","11_TxEnh3'","12_TxEnhW","13_EnhA1","14_EnhA2","15_EnhAF","16_EnhW1",
                            "17_EnhW2","18_EnhAc","19_DNase","20_ZNF/Rpts","21_Het","22_PromP","23_PromBiv","24_ReprPC","25_Quies"))

ggplot(imputed,aes(x=State,y=Count,fill=Imputed.State)) + geom_bar(position="fill",stat="identity") + ylab("Proportion") +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) + scale_fill_manual(values=imputed_states,name="State\n25-state model")
