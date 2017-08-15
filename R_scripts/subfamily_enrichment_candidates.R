# Subfamily enrichment candidate subfamiles
# See 7/11/2016, 8/29/2016, 9/29/2016, 11/22/2016, 11/27/2016, 5/23/2017, 5/24/2017, 5/29/2017, 5/30/2017, 6/15/2017

# Enhancer state enrichments to plot
enriched_7Enh_sigKW = intersect(enriched_subfamily$`7_Enh`,subfamily_state_kruskal_group[which(subfamily_state_kruskal_group$State == "7_Enh" & subfamily_state_kruskal_group$p_adjust < 0.05),]$Subfamily)

# Mesenchymal/epithelial/ENCODE/IMR90 sample grouping
c('E017', 'E023', 'E025', 'E026', 'E027', 'E028', 'E049', 'E052', 'E055', 'E056', 'E057', 'E058', 'E059', 'E061', 'E114', 'E117', 'E119', 'E120', 'E121', 'E122', 'E125', 'E126', 'E127', 'E128', 'E129')

# Subfamilies enriched in above sample group
epi_mesen_7Enh = c("LTR18A","LTR10D","LTR8","MER61C","MER61D","MER61E","HERVL18-int","LTR25","Charlie9","MER44C","MER44B","Tigger7","MER44D","MER121","LTR3B_","LTR26")

# 7_Enh subfamily-clustered samples
#sample_blood: E029,E030,E031,E032,E033,E034,E035,E036,E037,E038,E039,E040,E041,E042,E043,E044,E045,E046,E047,E048,E050,E051,E062,E093,E112,E116,E124,E115
#sample_brain:  E067,E068,E069,E070,E071,E072,E073,E074,E081,E082,E087
#sample_digestive:  E075,E077,E101,E102,E106,E109,E110,E066
#sample_ESCa:  E001,E002,E003,E004,E008,E011,E014,E015,E016,E018,E019,E020,E021,E024
#sample_ESCb:  E007,E009,E010,E012,E005,E013,E080,E084,E085,E086,E091,E099,E118,E123
#sample_mesen:  E017,E023,E025,E026,E049,E052,E055,E056,E059,E061,E114,E117,E120,E121,E122,E125,E126,E128,E129,E028,E057,E058,E119,E127
#sample_mesen_sub:  E028,E057,E058,E119,E127
#sample_other:  E006,E022,E027,E053,E054,E063,E065,E076,E078,E079,E083,E088,E089,E090,E092,E094,E095,E096,E097,E098,E100,E103,E104,E105,E107,E108,E111,E113
#sample_other_sub:  E006,E022,E079,E089,E094,E098

mer121_homer = as.data.frame(matrix(c("NF1-halfsite","Olig2","Atf3","AP-1","Fra1","BATF","BMAL1","Ets1-distal","ELF1","ETS",107,93,90,83,78,77,61,60,55,55),nrow=10,ncol=2))
colnames(mer121_homer) = c("TF","pvalue")
mer121_homer$pvalue = as.numeric(as.character(mer121_homer$pvalue))
mer121_GO = matrix(c("organ morphogenesis","embryonic morphogenesis","tissue morphogenesis","embryo development","tissue development","skeletal system development",48.9,39.7,38.4,37.8,35.9,35.5),nrow=6,ncol=2)
mer121_GO = as.data.frame(mer121_GO)
colnames(mer121_GO) = c("GO_Biological_Process","pvalue")
mer125_homer$pvalue = as.numeric(as.character(mer125_homer$pvalue))
colnames(mer125_homer) = c("TF","pvalue")
mer125_homer = as.data.frame(matrix(c("Pdx1","Nkx6.1","Lhx2","Lhx3","Lhx1",19,17,10,9,9),nrow=5,ncol=2))
x1_line_homer$pvalue = as.numeric(as.character(x1_line_homer$pvalue))
colnames(x1_line_homer) = c("TF","pvalue")
x1_line_homer = as.data.frame(matrix(c("EBF1","EBF","NeuroD1","Atoh1","Ascl1",16,11,9,8,6),nrow=5,ncol=2))

# Candidate subfamilies by state
candidate_1TssA = c('HERV1_I-int','HERV1_LTRa','HERV15-int','HERV3-int','HERVH-int','HERVKC4-int','HERVS71-int','LTR10B1','LTR10D','LTR12E','LTR12F','LTR13','LTR13_','LTR14C','LTR19C','LTR21A','LTR22A','LTR25','LTR26E','LTR2B','LTR2C','LTR30','LTR4','LTR43B','LTR46-int','LTR5_Hs','LTR6A','LTR7','LTR73','LTR75_1','LTR7Y','LTR9B','MamGyp-int','MER41G','MER50C','MER51A','MER51C','MER51D','MER57E3','MER61F')
candidate_2TssAFlnk = c('AmnSINE1','Charlie11','HERV15-int','HERV3-int','HERVIP10FH-int','HERVK11D-int','HERVK13-int','HERVKC4-int','HERVS71-int','LTR1','LTR10B','LTR10B1','LTR10D','LTR13_','LTR14','LTR14C','LTR19C','LTR21A','LTR22A','LTR23','LTR24','LTR24B','LTR25','LTR26B','LTR26E','LTR29','LTR2B','LTR2C','LTR30','LTR35A','LTR35B','LTR36','LTR3A','LTR3B','LTR4','LTR43B','LTR43-int','LTR44','LTR46','LTR48','LTR48B','LTR49','LTR5_Hs','LTR5B','LTR61','LTR62','LTR65','LTR6A','LTR7','LTR75_1','LTR77','LTR8','LTR8A','LTR9','LTR9B','MamGyp-int','MamRep1879','MamRep488','MER126','MER34C_','MER34C2','MER34D','MER41G','MER4A','MER4A1','MER4A1_','MER4C','MER4D','MER4D0','MER4D1','MER50C','MER51A','MER51B','MER51C','MER51D','MER57E3','MER61B','MER61F','MER67A','MER67B','MER67C','MER83A-int','MER84','MER88','MER92C','MER94B','PRIMA4_LTR','PrimLTR79','UCON26','UCON27','UCON28a','UCON29','UCON4')
candidate_3TxFlnk = c('HERV3-int','HERVH-int','HERVI-int','HERVIP10F-int','LFSINE_Vert','LTR10B1','LTR13','LTR7','MamGyp-int','MER21B','MER31-int','MER51C','MER67B','PRIMA41-int')
candidate_6EnhG = c('Charlie10a','Charlie10b','HERV4_I-int','HUERS-P3b-int','Kanga1b','L1M2a','L1P3b','LTR14','LTR19C','LTR34','LTR35A','LTR57','LTR5B','LTR65','LTR7','LTR70','LTR71B','LTR76','LTR86C','LTR87','LTR9B','MER105','MER107','MER21B','MER21-int','MER51B','MER52A','MER52D','MER67B','MER68','MER68B','MER6B','MLT1E1A-int','MLT1E2-int','Ricksha_b','Ricksha_c','Tigger11a','X7C_LINE')
candidate_7Enh = c('AmnSINE1','AmnSINE2','Charlie13a','Charlie13b','Eulor12','Eulor2A','Eulor5A','Eulor6D','Eulor8','Eulor9A','Eulor9C','HERV1_I-int','HERV1_LTRa','HERV1_LTRc','HERV1_LTRd','HERV1_LTRe','HERV15-int','HERV3-int','HERVK3-int','HERVL32-int','HERVS71-int','L1P4c','LFSINE_Vert','LTR10C','LTR10D','LTR13_','LTR13A','LTR14','LTR15','LTR18A','LTR18B','LTR19B','LTR19C','LTR2','LTR21A','LTR21B','LTR22','LTR22A','LTR22B','LTR23','LTR24','LTR24B','LTR24C','LTR26B','LTR26E','LTR29','LTR2B','LTR2C','LTR30','LTR34','LTR35','LTR35A','LTR35B','LTR36','LTR38C','LTR39','LTR3A','LTR3B','LTR3B_','LTR43','LTR43-int','LTR44','LTR46','LTR48','LTR48B','LTR49','LTR5_Hs','LTR51','LTR59','LTR5B','LTR61','LTR64','LTR65','LTR6A','LTR6B','LTR7','LTR71B','LTR72','LTR75_1','LTR76','LTR77','LTR7C','LTR8','LTR9','LTR9B','MamSINE1','MER121','MER123','MER126','MER129','MER130','MER131','MER134','MER21A','MER31A','MER34A1','MER34B','MER34C','MER34C_','MER34C2','MER34D','MER39B','MER41A','MER41B','MER41C','MER41G','MER44B','MER44C','MER44D','MER48','MER4A','MER4A1','MER4A1_','MER4B','MER4C','MER4D','MER4D0','MER4D1','MER4E','MER4E1','MER50C','MER51B','MER51C','MER51D','MER57E3','MER61B','MER61C','MER61D','MER61E','MER61F','MER67B','MER67C','MER67D','MER72B','MER74A','MER74B','MER83','MER83A-int','MER91A','MER9B','Merlin1_HS','PABL_A','PABL_B','PRIMA4_LTR','PrimLTR79','Tigger1a_Mars','UCON10','UCON11','UCON12','UCON13','UCON14','UCON17','UCON2','UCON20','UCON26','UCON27','UCON29','UCON31','UCON4','UCON6','UCON8','UCON9','X1_LINE','X2_LINE','X3_LINE','X9_LINE')

# Statistics for candidate subfamilies
test = merge(merge(aggregate(data=subfamily_state_sample_filter[which(subfamily_state_sample_filter$Subfamily %in% candidate_1TssA & subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "1_TssA"),],Members~Subfamily,min),aggregate(data=subfamily_state_sample_filter[which(subfamily_state_sample_filter$Subfamily %in% candidate_1TssA & subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "1_TssA"),],Members~Subfamily,max),by=c("Subfamily")),aggregate(data=subfamily_state_sample_filter[which(subfamily_state_sample_filter$Subfamily %in% candidate_1TssA & subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "1_TssA"),],Members~Subfamily,median),by=c("Subfamily"))
colnames(test) = c("Subfamily","Min","Max","Median")
test = merge(test,aggregate(data=potential_TEother_state[which(potential_TEother_state$subfamily %in% candidate_1TssA),],X1_TssA~subfamily,function(x) sum(x > 0)),by.x=c("Subfamily"),by.y=c("subfamily"))
test = merge(test,subfamily_state_sample_counts[which(subfamily_state_sample_counts$State == "1_TssA"),c(3,5)],by=c("Subfamily"))
colnames(test)[5:6] = c("Members_ever","Sample_enriched")
test = merge(test,rmsk_TEother_stats_subfamily[,3:4],by=c("Subfamily"))
colnames(test)[7] = "Members"
test = merge(test,aggregate(data=potential_TEother_state[which(potential_TEother_state$subfamily %in% candidate_1TssA),],X1_TssA~subfamily,function(x) sum(x > 5)),by.x=c("Subfamily"),by.y=c("subfamily"))

# Statistics for candidate subfamilies
test = merge(merge(aggregate(data=subfamily_state_sample_filter[which(subfamily_state_sample_filter$Subfamily %in% candidate_7Enh & subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "7_Enh"),],Members~Subfamily,min),aggregate(data=subfamily_state_sample_filter[which(subfamily_state_sample_filter$Subfamily %in% candidate_7Enh & subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "7_Enh"),],Members~Subfamily,max),by=c("Subfamily")),aggregate(data=subfamily_state_sample_filter[which(subfamily_state_sample_filter$Subfamily %in% candidate_7Enh & subfamily_state_sample_filter$Enrichment > 1.5 & subfamily_state_sample_filter$State == "7_Enh"),],Members~Subfamily,median),by=c("Subfamily"))
colnames(test) = c("Subfamily","Min","Max","Median")
test = merge(test,aggregate(data=potential_TEother_state[which(potential_TEother_state$subfamily %in% candidate_7Enh),],X7_Enh~subfamily,function(x) sum(x > 0)),by.x=c("Subfamily"),by.y=c("subfamily"))
test = merge(test,subfamily_state_sample_counts[which(subfamily_state_sample_counts$State == "7_Enh"),c(3,5)],by=c("Subfamily"))
colnames(test)[5:6] = c("Members_ever","Sample_enriched")
test = merge(test,rmsk_TEother_stats_subfamily[,3:4],by=c("Subfamily"))
colnames(test)[7] = "Members"
test = merge(test,aggregate(data=potential_TEother_state[which(potential_TEother_state$subfamily %in% candidate_7Enh),],X7_Enh~subfamily,function(x) sum(x > 5)),by.x=c("Subfamily"),by.y=c("subfamily"))

# Individual TEs in candidate subfamilies in state by sample
candidate_1TssA_indv = read.table("candidate_1TssA_indv.txt",sep='\t')
candidate_1TssA_indv = unique(candidate_1TssA_indv)
colnames(candidate_1TssA_indv) = c("chromosome","start","stop","subfamily","class","family","strand","State","Overlap","Sample")

candidate_7Enh_indv = read.table("enrichment/candidate_7Enh_indv.txt",sep='\t')
colnames(candidate_7Enh_indv) = colnames(candidate_1TssA_indv)
