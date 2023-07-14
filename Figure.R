rm(list=ls())
options( stringsAsFactors = F )
setwd("H:\\BMC_top\\Figure_1")
#-----------------------------------------------------------------
library(maftools)
library(dplyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(impute)
library(wateRmelon)
library(dplyr)
library(minfi)
library(ComplexHeatmap)
library(matrixStats)
library(circlize)
library(RColorBrewer)
library(GetoptLong)
library(GenomicRanges)
library(ConsensusClusterPlus)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(reshape2)
library(tableone)
require(ggplot2)
require(ggpubr)
library(tableone)
library(ggtern)
#---------------------------------------------------------------
#获取已知低表达的结肠癌抑癌基因
#Min Zhao, Jingchun Sun, Zhongming Zhao. 
#TSGene: a web resource for tumor suppressor genes.
#Nucleic Acids Research, 2013 Jan;41(Database issue):D970-6.
TSG_colon <- colon.TSGs$Gene
length(TSG_colon)
#总共535个TSG
#----------------------------------------------------------------
#获取450k芯片注释文件
#官方注释发表文章:Validation of a DNA methylation microarray 
#for 450,000 CpG sites in the human genome
#将空白位置转换为intergenic
ann450k2$island_location[ann450k2$island_location == ""] <- "Intergenic"
table(ann450k2$island_location)
#1stExon      3'UTR      5'UTR      Body   Intergenic    TSS1500     TSS200 
# 10810      15379      49521     150147     119652      77378      62625 
table(ann450k2$Relation_to_Island)
#Island N_Shelf N_Shore OpenSea S_Shelf S_Shore 
#150254   24844   62870  176047   22300   49197
#获取注释完毕,与官方注释一致
#---------------------------------------------------------------------------
#获取TSG的甲基化位点
TSG_ann450k <- ann450k2 %>%
  dplyr::filter(gene_id %in% TSG_colon)
#只有528个gene对应到甲基化位点
length(unique(TSG_ann450k$gene_id))

#获取TSG promoter区域的甲基化位点
TSG_ann450k_promoter <-TSG_ann450k[TSG_ann450k$island_location == "TSS200"|
                                     TSG_ann450k$island_location =="TSS1500",] 
#获取4697个probes
dim(TSG_ann450k_promoter)
table(TSG_ann450k_promoter$island_location)
#TSS1500  TSS200 
#2832     1865
#---------------------------------------------------------------------------
#导入结肠癌甲基化矩阵
#获取甲基化探针
row.names(Data) <- Des[,1] 
#526732 probes  353 samples
dim(Data)  
#预处理甲基化芯片
#第一步:删除缺失值，剩余409441 probes
Data1 <- na.omit(Data)
dim(Data1)
#第二步:删除重复的probes，有34636个重复探针，留下最大值，剩余374805
sum(duplicated(rownames(Data1)))
Data1 <- Data1[order(rowMeans(Data1)),]
matData1 = Data1[!duplicated(row.names(Data1)),]
dim(matData1)
#第三步:去除性染色体上的9029个探针,剩余365776个探针
#删除探针落在X Y 染色体上
keep <- !(row.names(matData1) %in% ann450k2$Name[ann450k2$chr %in% 
                                                  c("chrX","chrY")])
table(keep)
matData1 <- matData1[keep,]
#第四步:remove probes with SNPs at CpG site,剩余364236个探针
#将beta矩阵转换为minfi包的GenomicRatioSet格式
grset=makeGenomicRatioSetFromMatrix(matData1,what="Beta")
class(grset)
mSetSqFlt <- dropLociWithSnps(grset)
dim(mSetSqFlt)
#第五步:exclude cross reactive probes,剩余349940个探针
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)
mSetSqFlt <- mSetSqFlt[keep,] 
dim(mSetSqFlt)
#预处理完毕,将GenomicRatioSet格式转换为被他矩阵
bVals <- getBeta(mSetSqFlt)
dim(bVals)
#数据质控完毕，获得beta矩阵：bVals，349940个probes，353个样本
#---------------------------------------------------------------------------
#获取临床表达矩阵
metadata <- data.frame(colnames(bVals))

for (i in 1:length(metadata[,1])) {
  
  num <- as.numeric(substring(metadata[i,1],14,15))
  
  if (num %in% seq(1,9)) {metadata[i,2] <- "Tumor"}
  
  if (num %in% seq(10,29)) {metadata[i,2] <- "Normal"}
  
}
names(metadata) <- c("TCGA_id","sample")
metadata$sample <- as.factor(metadata$sample)
metadata %>% dplyr::group_by(sample) %>% summarise(n())
#癌旁 38个，癌症 315个
#--------------------------------------------------------------------------
#正常样本beta矩阵
#检查barcode是否对齐
sum(metadata$TCGA_id == colnames(bVals))
Normal_barcode <- colnames(bVals)[metadata$sample =="Normal"]
bVals_df <- as.data.frame(bVals)
#获取正常矩阵:Normal_bVals,349940个probes,38个样本
Normal_bVals <- bVals_df %>%
  dplyr::select(Normal_barcode)
dim(Normal_bVals)
#肿瘤样本矩阵:Tumor_bVals,349940个probes,315个样本
Tumor_barcode <- colnames(bVals)[metadata$sample =="Tumor"]
Tumor_bVals <- bVals_df %>%
  dplyr::select(Tumor_barcode)
dim(Tumor_bVals)

#减少肿瘤样本中的重复
colnames(Tumor_bVals) <- substr(colnames(Tumor_bVals),1,15)
#去除重复值后，剩余298个样本
Tumor_bVals = Tumor_bVals[,!duplicated(colnames(Tumor_bVals))]
dim(Tumor_bVals)
Tumor_bVals <- Tumor_bVals[,substr(colnames(Tumor_bVals),14,15) == "01"]
#删除不是01结尾的样本，剩余296个样本
dim(Tumor_bVals)
#获得肿瘤矩阵：Tumor_bVals，349940个探针，296个样本
#--------------------------------------------------------------------------
class(Tumor_bVals)
Tumor_bVals$probe <- row.names(Tumor_bVals)
Tumor_bVals_promoter <- Tumor_bVals %>%
  filter(probe %in% TSG_ann450k_promoter$Name)
#从4697个probes减少到3783
dim(Tumor_bVals_promoter)

#肿瘤组的promoter矩阵
row.names(Tumor_bVals_promoter) <- Tumor_bVals_promoter$probe
Tumor_bVals_promoter <- Tumor_bVals_promoter[,-297]
dim(Tumor_bVals_promoter)
#标准差大于0.2,获得481个probes：final_Tumor
Tumor_bVals_promoter$SD <- apply(Tumor_bVals_promoter,1,sd)
sum(Tumor_bVals_promoter$SD>0.2)
inclusion_Tumor <- (Tumor_bVals_promoter$SD>0.2)
final_Tumor <- row.names(Tumor_bVals_promoter)[inclusion_Tumor]
length(final_Tumor)

#正常组的平均值小于0.05，获得1493probes：final_Normal
class(Normal_bVals)
Normal_bVals$probe <- row.names(Normal_bVals)
Normal_bVals_promoter <- Normal_bVals %>%
  filter(probe %in% TSG_ann450k_promoter$Name)
dim(Normal_bVals_promoter)
row.names(Normal_bVals_promoter) <- Normal_bVals_promoter$probe
Normal_bVals_promoter <- Normal_bVals_promoter[,-39]
Normal_bVals_promoter$average <- apply(Normal_bVals_promoter,1,mean)
inclusion_Normal <- (Normal_bVals_promoter$average < 0.05)
final_Normal <- row.names(Normal_bVals_promoter)[inclusion_Normal ]
length(final_Normal )
#两者取交集
clust_probe <- intersect(final_Tumor,final_Normal)
#获得171个聚类probes
length(clust_probe)
#获取171个probes的完整注释
clust_ann450k <- ann450k2 %>%
  dplyr::filter(Name %in% clust_probe)
dim(clust_ann450k)
write.table(clust_ann450k,"clust_probe.csv",sep=",")
table(clust_ann450k$island_location)
#TSS1500  TSS200 
#  88      83 
table(clust_ann450k$Relation_to_Island)
#Island N_Shore OpenSea S_Shore 
# 142      13       7       9
#-----------------------------------------------------------------------
#获取聚类的肿瘤样本
Tumor_clust <- Tumor_bVals %>%
  filter(probe %in% clust_probe)
row.names(Tumor_clust) <- Tumor_clust$probe
#获得肿瘤聚类矩阵Tumor_clust
Tumor_clust <- Tumor_clust[,-297]
#171个probes,296个样本
dim(Tumor_clust)

###构建癌旁聚类矩阵
Normal_clust <- Normal_bVals %>%
  filter(probe %in% clust_probe)
row.names(Normal_clust) <- Normal_clust$probe
Normal_clust <- Normal_clust[,-39]
#获取171个probes,38个样本
dim(Normal_clust)
#讲数据框转换为matrix
Tumor_clust1 <- as.matrix(Tumor_clust)
results<-ConsensusClusterPlus(Tumor_clust1,maxK=10,pItem=0.8,
                              clusterAlg = 'km',distance='euclidean',
                              plot="pdf",seed= 200000)
#k值取3
clust_result <- results[[3]][["consensusClass"]]
table(clust_result)
#  1   2   3 
# 99 141  56 
56/296
#对结果进行排序
clust_result <- clust_result[order(clust_result)]
clust_result <- as.data.frame(clust_result)
#将肿瘤样本顺序按照聚类结果排序
Tumor_clust1 <- Tumor_clust %>%
  dplyr::select(row.names(clust_result))
#查看结果是否对齐
colnames(Tumor_clust1) ==row.names(clust_result )
#将聚类结果转换为字符串
clust_result$clust_result[clust_result$clust_result ==1] <- "I"
clust_result$clust_result[clust_result$clust_result ==2] <- "II"
clust_result$clust_result[clust_result$clust_result ==3] <- "III" 
#----------------------------------------------------------------
#提取聚类probes注释子集
clust_ann450k_final <- clust_ann450k[,c(4,19,34)]
#----------------------------------------------------------------
#导入临床信息
clinical <- COAD_clinicalMatrix
#获取296个病人的临床信息
clinical <- clinical %>%
  filter(sampleID %in% colnames(Tumor_clust1))
dim(clinical)

#更改肿瘤位置的临床信息
table(clinical$anatomic_neoplasm_subdivision)
clinical$anatomic_neoplasm_subdivision[clinical$anatomic_neoplasm_subdivision =="Ascending Colon"] <- "Right"
clinical$anatomic_neoplasm_subdivision[clinical$anatomic_neoplasm_subdivision =="Cecum"] <- "Right"
clinical$anatomic_neoplasm_subdivision[clinical$anatomic_neoplasm_subdivision =="Hepatic Flexure"] <- "Right"
clinical$anatomic_neoplasm_subdivision[clinical$anatomic_neoplasm_subdivision =="Transverse Colon"] <- "Right"
clinical$anatomic_neoplasm_subdivision[clinical$anatomic_neoplasm_subdivision =="Splenic Flexure"] <- "Left"
clinical$anatomic_neoplasm_subdivision[clinical$anatomic_neoplasm_subdivision =="Sigmoid Colon"] <- "Left"
clinical$anatomic_neoplasm_subdivision[clinical$anatomic_neoplasm_subdivision =="Rectosigmoid Junction"] <- "Left"
clinical$anatomic_neoplasm_subdivision[clinical$anatomic_neoplasm_subdivision =="Descending Colon"] <- "Left"
clinical$anatomic_neoplasm_subdivision[clinical$anatomic_neoplasm_subdivision ==""] <- NA
clinical$anatomic_neoplasm_subdivision[clinical$anatomic_neoplasm_subdivision =="[Discrepancy]"] <- NA
table(clinical$anatomic_neoplasm_subdivision)
#Left Right 
#103   176

#更改生存状态的信息
table(clinical$vital_status)
clinical$vital_status[clinical$vital_status == "LIVING"] <- "Alive"
clinical$vital_status[clinical$vital_status == "DECEASED"] <- "Dead"
clinical$vital_status[clinical$vital_status ==""] <- NA
table(clinical$vital_status)
#Alive  Dead 
#224    70


#更改性别信息
table(clinical$gender)
clinical$gender[clinical$gender =="MALE"] <- "Male"
clinical$gender[clinical$gender =="FEMALE"] <- "Female"
clinical$gender[clinical$gender ==""] <- NA
table(clinical$gender)
#Female   Male 
#136    158 

#更改结肠息肉病史信息
table(clinical$history_of_colon_polyps)
clinical$history_of_colon_polyps[clinical$history_of_colon_polyps =="YES"] <- "Yes"
clinical$history_of_colon_polyps[clinical$history_of_colon_polyps =="NO"] <- "No"
clinical$history_of_colon_polyps[clinical$history_of_colon_polyps ==""] <- NA
table(clinical$history_of_colon_polyps)
#No Yes 
#172  53 

#更改pathologic_M信息
table(clinical$pathologic_M)
clinical$pathologic_M[clinical$pathologic_M=="M1a"] <- "M1"
clinical$pathologic_M[clinical$pathologic_M=="M1b"] <- "M1"
clinical$pathologic_M[clinical$pathologic_M==""] <- NA
table(clinical$pathologic_M)
# M0  M1  MX 
#198  41  50 

#更改pathologic_T信息
table(clinical$pathologic_T)
clinical$pathologic_T[clinical$pathologic_T == "T4a"] <- "T4"
clinical$pathologic_T[clinical$pathologic_T == "T4b"] <- "T4"
clinical$pathologic_T[clinical$pathologic_T == "Tis"] <- NA
clinical$pathologic_T[clinical$pathologic_T == ""] <- NA
table(clinical$pathologic_T)
#T1  T2  T3  T4 
#7  43 203  40 

#更改pathologic_N信息
table(clinical$pathologic_N)
clinical$pathologic_N[clinical$pathologic_N == "N1a"] <- "N1"
clinical$pathologic_N[clinical$pathologic_N == "N1b"] <- "N1"
clinical$pathologic_N[clinical$pathologic_N == "N1c"] <- "N1"
clinical$pathologic_N[clinical$pathologic_N == "N2a"] <- "N2"
clinical$pathologic_N[clinical$pathologic_N == "N2b"] <- "N2"
clinical$pathologic_N[clinical$pathologic_N == ""] <- NA
table(clinical$pathologic_N)
#N0  N1  N2 
#171  74  49

#更改pathologic_stage信息
table(clinical$pathologic_stage)
clinical$pathologic_stage[clinical$pathologic_stage =="Stage IA"] <- "Stage I"
clinical$pathologic_stage[clinical$pathologic_stage =="Stage IIA"] <- "Stage II"
clinical$pathologic_stage[clinical$pathologic_stage =="Stage IIB"] <- "Stage II"
clinical$pathologic_stage[clinical$pathologic_stage =="Stage IIC"] <- "Stage II"
clinical$pathologic_stage[clinical$pathologic_stage =="Stage IIIA"] <- "Stage III"
clinical$pathologic_stage[clinical$pathologic_stage =="Stage IIIB"] <- "Stage III"
clinical$pathologic_stage[clinical$pathologic_stage =="Stage IIIC"] <- "Stage III"
clinical$pathologic_stage[clinical$pathologic_stage =="Stage IVA"] <- "Stage IV"
clinical$pathologic_stage[clinical$pathologic_stage =="Stage IVB"] <- "Stage IV"
clinical$pathologic_stage[clinical$pathologic_stage =="[Discrepancy]"] <- NA
clinical$pathologic_stage[clinical$pathologic_stage ==""] <- NA
table(clinical$pathologic_stage)
#Stage I  Stage II Stage III  Stage IV 
#44       115        85        41 

#更改venous_invasion信息
table(clinical$venous_invasion)
clinical$venous_invasion[clinical$venous_invasion == "YES"] <- "Yes"
clinical$venous_invasion[clinical$venous_invasion == "NO"] <- "No"
clinical$venous_invasion[clinical$venous_invasion == ""] <- NA
table(clinical$venous_invasion)
#No Yes 
#196  60

#更改名字
clinical1 <- clinical %>%
  dplyr::rename(PatientID = sampleID,
                Age = age_at_initial_pathologic_diagnosis,
                Anatomic_location = anatomic_neoplasm_subdivision)

#数据框构建完毕，clinical1，clust_ann450k_final，Tumor_clust1，clust_result,Normal_clust
#第一步,过关
row.names(Tumor_clust1) == row.names(Normal_clust)
#第二步,过关
colnames(Tumor_clust1) == row.names(clust_result)
#第三步，clinical4过关
clinical1$PatientID == colnames(Tumor_clust1)
clinical2 <- t(clinical1)
colnames(clinical2) <- clinical2[1,]
clinical3 <- as.data.frame(clinical2)
clinical4 <- clinical3 %>%
  dplyr::select(colnames(Tumor_clust1))
clinical4 <- t(clinical4)
clinical4 <- as.data.frame(clinical4)
clinical4$PatientID == colnames(Tumor_clust1)
#第四步clust_ann450k_final1过关
clust_ann450k_final$Name == rownames(Tumor_clust1)
clust_ann450k_final <- t(clust_ann450k_final)
colnames(clust_ann450k_final) <- clust_ann450k_final[1,]
clust_ann450k_final <- as.data.frame(clust_ann450k_final)
clust_ann450k_final1 <- clust_ann450k_final %>%
  dplyr::select(rownames(Tumor_clust1))
clust_ann450k_final1 <- t(clust_ann450k_final1)
clust_ann450k_final1 <- as.data.frame(clust_ann450k_final1)
clust_ann450k_final1$Name == rownames(Tumor_clust1)

#更改临床数据中的年龄
clinical4$Age <- as.numeric(clinical4$Age)
sum(clinical4$Age>65)
clinical4$Age <- ifelse(clinical4$Age >65,"High","Low")
table(clinical4$Age)
#High  Low 
#153  141 
#检查完毕，clinical4，clust_ann450k_final1，Tumor_clust1，clust_result,Normal_clust
#------------------------------------------------------------------
######绘制热图
ha = HeatmapAnnotation(Cluster = clust_result$clust_result,
                       Stage = clinical4$pathologic_stage,
                       Gender = clinical4$gender,
                       Anatomic_location = clinical4$Anatomic_location,
                       vital_status = clinical4$vital_status,
                       Age = clinical4$Age,
                       History_of_polyps = clinical4$history_of_colon_polyps,
                       Venous_invasion = clinical4$venous_invasion,
                       col = list(Cluster = structure(names = c("I", "II", "III", "IV"), brewer.pal(4, "Set1")),
                                  Stage = structure(names = c("Stage I", "Stage II", "Stage III", "Stage IV"), 
                                                    brewer.pal(4, "Set1")),
                                  Gender = structure(names = c("Female", "Male"), c("#377EB8","#4DAF4A")),
                                  Anatomic_location = structure(names = c("Right", "Left"),c("#377EB8","#4DAF4A")),
                                  vital_status = structure(names = c("Alive", "Dead"),c("white","#377EB8")),
                                  Age = structure(names = c("High", "Low"), c("#E41A1C","white")),
                                  History_of_polyps = structure(names = c("Yes", "No"),c("#377EB8", "white")),
                                  Venous_invasion = structure(names = c("Yes", "No"),c("#E41A1C","white"))),
                       na_col = "grey",
                       show_legend = c(TRUE, TRUE, TRUE,TRUE,F,F,F,F),
                       annotation_height = unit(c( 5, 5, 5,5,5,5, 5, 5), "mm"),
                       annotation_legend_param = list(
                         Cluster = list(title = "Cluster"),
                         Stage = list(title = "Stage"),
                         Anatomic_location = list(title = "Lesion_location"),
                         Gender  = list(title = "Gender")))

col_fun = colorRamp2(c(0, 0.5, 1), c("#377EB8", "white", "#E41A1C"))


ht_list = Heatmap(Tumor_clust1, col = col_fun, name = "Methylation",
                  clustering_distance_rows = "euclidean", row_dend_reorder = TRUE,
                  cluster_columns = F, cluster_rows = F,
                  show_row_dend = FALSE, show_column_dend = FALSE,
                  show_row_names = FALSE, show_column_names = FALSE,
                  bottom_annotation = ha, column_title = qq("COAD samples (n = @{ncol(Tumor_clust1)})"),
                  column_title_gp = gpar(font = 2),
                  row_title_gp = gpar(col = "#FFFFFF00"))+ 
  Heatmap(Normal_clust, col = col_fun, cluster_rows = F,show_row_names = FALSE, show_column_names = FALSE, 
          show_column_dend = FALSE, column_title = "Normal",
          column_title_gp = gpar(font = 2),
          show_heatmap_legend = FALSE, width = unit(1, "cm"))+
  Heatmap(clust_ann450k_final1$Relation_to_Island,
          name = "Relation to CGI", cluster_rows = F,
          show_row_names = FALSE, show_column_names = FALSE, width = unit(5, "mm"),
          col = structure(names = c("Island", "N_Shore", "OpenSea", "S_Shelf","S_Shore"), 
                          c("#377EB8","#4DAF4A","#984EA3","#FF7F00","#E41A1C")))+
  Heatmap(clust_ann450k_final1$island_location,
          name = "Relation to TSS", cluster_rows = F,
          show_row_names = FALSE, show_column_names = FALSE, 
          width = unit(5, "mm"),
          col = structure(names = c("TSS200", "TSS1500"), 
                          c("#377EB8","#4DAF4A")))
draw(ht_list, annotation_legend_side = "left", heatmap_legend_side = "left")

annotation_titles = c(Cluster = "Cluster",
                      Stage = "Pathologic stage",
                      Gender = "Gender",
                      Anatomic_location = "Anatomic_location",
                      vital_status = "Staus_dead",
                      Age = "Age_high",
                      History_of_polyps = "History of polyps",
                      Venous_invasion = "Venous invasion")
for(an in names(annotation_titles)) {
  decorate_annotation(an, {
    grid.text(annotation_titles[an], unit(-2, "mm"), just = "right",gp = gpar(font = 2))
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}

decorate_heatmap_body("Relation to CGI", slice = NULL, {
  grid.text("Relation to CGI", y = unit(-2, "mm"), rot = 90, just = "right",gp = gpar(font = 2))
})
decorate_heatmap_body("Relation to TSS", slice = NULL, {
  grid.text("Relation to TSS", y = unit(-1.5, "mm"), 
            rot = 90, just = "right",gp = gpar(font = 2))
})
#heatmap绘制完毕，保存图像14:8
#--------------------------------------------------------------------------
#突变数据,按照maftools包的要求将名字更改
colnames(clinical4)[1] <- "Tumor_Sample_Barcode"
clinical5 <- clinical4 %>%
  dplyr::rename(Overall_Survival_Status = OS, time = OS.time)

#合并聚类和临床矩阵
class(clust_result)
clust_result$Tumor_Sample_Barcode <- row.names(clust_result)

clinical6 <- inner_join(clinical5,clust_result , by= "Tumor_Sample_Barcode")
clinical6 <- clinical6 %>%
  dplyr::rename(Cluster = clust_result)
#按照maftools包的要求，更改样本的名称
clinical6$Tumor_Sample_Barcode <- substr(clinical6$Tumor_Sample_Barcode,1,12)
#导入下载好的突变数据
load("H:/Genome_Medicine/COADmutation.RData")
maf <- read.maf(maf = COADmut, clinicalData = clinical6, isTCGA = T)

#获取与临床样本一致的子集
maf_clinical <- subsetMaf(maf = maf, tsb = clinical6$Tumor_Sample_Barcode, 
                          mafObj = TRUE)

#只有291个样本有突变数据


#绘制oncoplot图，保存图像14:8
col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
fabcolors = RColorBrewer::brewer.pal(n = 3,name = 'Set1')

names(fabcolors) = c("I", "II", "III")
fabcolors = list(Cluster = fabcolors)

oncoplot(maf = maf_clinical,colors = col,
         top = 60,
         legendFontSize = 8,
         clinicalFeatures = 'Cluster',
         sortByAnnotation = TRUE, 
         annotationColor = fabcolors)
decorate_annotation("Cluster", {
  grid.text("Cluster", x = unit(-0.5, "mm"), just = "right")
})

save("maf","maf_clinical","clinical6",file = "mutation.RData")
#---------------------------------------------------------------------
#Go on
setwd("H:\\BMC_top\\Figure_1")

table(maf@clinical.data$Cluster)
#I  II III 
#99 136  56

#g根据临床聚类结果取突变子集
ClusterI_ID <- clinical6$Tumor_Sample_Barcode[clinical6$Cluster == "I"]
maf_I <- subsetMaf(maf = maf, tsb = ClusterI_ID, mafObj = TRUE)

ClusterII_ID <- clinical6$Tumor_Sample_Barcode[clinical6$Cluster == "II"]
maf_II <- subsetMaf(maf = maf, tsb = ClusterII_ID, mafObj = TRUE)


ClusterIII_ID <- clinical6$Tumor_Sample_Barcode[clinical6$Cluster == "III"]
maf_III <- subsetMaf(maf = maf, tsb = ClusterIII_ID, mafObj = TRUE)
length(maf_II@clinical.data$Cluster)
##获得maf_clinical各个cluster的子集:ClusterI_ID,ClusterII_ID,ClusterIII_ID

#获取每个cluster子集中基因的突变频率
mutation_I <- as.data.frame(getGeneSummary(maf_I))
mutation_II <- as.data.frame(getGeneSummary(maf_II))
mutation_III <- as.data.frame(getGeneSummary(maf_III))

#找出各组突变top30，并取并集
top30 <- unique(c(mutation_I$Hugo_Symbol[1:30],mutation_II$Hugo_Symbol[1:30],
                  mutation_III$Hugo_Symbol[1:30]))
#共获得55个genes
length(top30)
#在各组突变频率数据中取出top30的子集

mutation_I_1 <- mutation_I %>%
  filter(Hugo_Symbol %in% top30)  
mutation_II_1 <- mutation_II %>%
  filter(Hugo_Symbol %in% top30)  
mutation_III_1 <- mutation_III %>%
  filter(Hugo_Symbol %in% top30) 
names(mutation_III_1)

mutation_III_1 <- mutation_III_1[,c(1,12)]
mutation_I_1 <- mutation_I_1[,c(1,12)]
mutation_II_1 <- mutation_II_1[,c(1,12)]

mutation_I_1$Cluster <- "I"
mutation_II_1$Cluster <- "II"
mutation_III_1$Cluster <- "III"
#获得整个top genes 的突变频率
mutation_TOTAL <- rbind(mutation_I_1,mutation_II_1,mutation_III_1)
#将数据变形
mutation_TOTAL <- spread(mutation_TOTAL,Cluster,MutatedSamples)

mutation_TOTAL$I_percent <- mutation_TOTAL$I/99
mutation_TOTAL$II_percent <- mutation_TOTAL$II/136
mutation_TOTAL$III_percent <- mutation_TOTAL$III/56

mutation_TOTAL$P <- NA
#颠换行与列
row.names(mutation_TOTAL) <- mutation_TOTAL$Hugo_Symbol
mutation_TOTAL <- mutation_TOTAL[,-1]
mutation_TOTAL1 <- t(mutation_TOTAL)

mutation_TOTAL1 <- as.data.frame(mutation_TOTAL1)
mutation_TOTAL1[c(4:6),] <- mutation_TOTAL1[c(4:6),]*100
#对数据进行卡方检验
x <- c()
for(i in 1:ncol(mutation_TOTAL1)){
  x <- chisq.test(mutation_TOTAL1[4:6,i])
  mutation_TOTAL1[7,i] <- x$p.value
}

mutation_TOTAL2 <- as.data.frame(t(mutation_TOTAL1))
write.table(mutation_TOTAL2,"mutation.csv",sep=",")

#----------------------------------------------------------------
###肿瘤突变负荷的计算
#讲maf结果转换为矩阵

#获取每个cluster中每个样本的突变总数
TMB_I <- as.data.frame(maf_I@variants.per.sample)
TMB_II <- as.data.frame(maf_II@variants.per.sample)
TMB_III <- as.data.frame(maf_III@variants.per.sample)

#人类有35*1Mbp个碱基，TMB是每1MB中突变的数目
TMB_I$TMB <- TMB_I$Variants%/%35
TMB_II$TMB <- TMB_II$Variants%/%35
TMB_III$TMB <- TMB_III$Variants%/%35

TMB_I$Cluster <- "I"
TMB_II$Cluster <- "II"
TMB_III$Cluster <- "III"
#将数据合并
TMB_total <- rbind(TMB_I,TMB_II,TMB_III)
#将类转换为因子
TMB_total$Cluster <- factor(TMB_total$Cluster)

my_comparisons <- list(c("I", "II"), c("I", "III"), 
                       c("II","III"))
result <- as.data.frame(compare_means(TMB~Cluster, data=TMB_total,method = "wilcox.test"))
#将结果保存
write.table(result,"TMB_CLUSTER.csv",sep=",")
write.table(TMB_total,"TMB_result.csv",sep=",")

#绘图
p <- ggplot(TMB_total, aes(x = Cluster, y = TMB, color = Cluster))+
  geom_boxplot(outlier.color = NA)+
  stat_compare_means(method = "kruskal.test", label.x = 0.56,label.y = 10) + #label的位置
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "II")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_blank()) + #去除外层边框
  theme(axis.line = element_line(colour = "black")) + #沿坐标轴显示直线
  xlab("Cluster") +
  theme(axis.title.x=element_text(size=12,face="bold"), 
        axis.title.y=element_text(size=12,face="bold"), 
        axis.text=element_text(size=10,face="bold"), 
        legend.text = element_text( size = 10, face = 'bold'))+
  guides(color=FALSE)


p +  ylim(0,60)+geom_dotplot(binaxis = "y", #沿y轴堆积，并沿着x轴分组
                             binwidth = 0.5, #最大组距
                             dotsize = 1.2, #点的大小
                             #如果点太多，两组叠在一起，就需要运行下面这行把它们分开
                             #stackgroups = T, binpositions="all",
                             stackdir = "center")



TMB_total$TMB_classific <- TMB_total$TMB
N1 <- TMB_total$TMB_classific>=20
N2 <- TMB_total$TMB_classific < 6
N3 <- TMB_total$TMB_classific >=6 & TMB_total$TMB_classific < 20
TMB_total$TMB_classific[N1] <- "High"
TMB_total$TMB_classific[N2] <- "Low"
TMB_total$TMB_classific[N3] <- "Mid"
table(TMB_total$TMB_classific,TMB_total$Cluster)
#按照指南分组的情况
#      I  II III
#High  12  12  31
#Low   84 118  24
#Mid    3   6   1
#-----------------------------------------------------------------------------
#将数据中的临床信息储存起来
clinical6$pathologic_T[clinical6$pathologic_T == "T1"] <- "T1+T2"
clinical6$pathologic_T[clinical6$pathologic_T == "T2"] <- "T1+T2"
table(clinical6$pathologic_T)

vars <- names(clinical6)[6:15]
strata <- "Cluster"
dim(clinical6)
dataFrame <- clinical6[,6:16]
tableOne <- CreateTableOne(vars = vars, strata = strata, data = dataFrame)
print(tableOne,quote = T,noSpaces = T)
#-----------------------------------------------------------------------------
#获取PI3K通路的突变比例
PI3K_gene <- c("IGF2","ERBB2","ERBB3","IRS2","PTEN",
               "NRAS","KRAS","PIK3R1","PIK3CA","BRAF")

mutation_III_PI3K <- mutation_III %>%
  filter(Hugo_Symbol %in% PI3K_gene)
mutation_II_PI3K <- mutation_II %>%
  filter(Hugo_Symbol %in% PI3K_gene)
mutation_I_PI3K <- mutation_I %>%
  filter(Hugo_Symbol %in% PI3K_gene)
#获得PI3K通路的突变数据
mutation_I_PI3K <- mutation_I_PI3K[,c(1,12)]
mutation_III_PI3K <- mutation_III_PI3K[,c(1,12)]
mutation_II_PI3K <- mutation_II_PI3K[,c(1,12)]

mutation_III_PI3K$Cluster <- "CIMP-H"
mutation_II_PI3K$Cluster <- "CIMP-L"
mutation_I_PI3K$Cluster <- "non-CIMP"
mutation_PI3K <- rbind(mutation_I_PI3K,mutation_III_PI3K,mutation_II_PI3K)

mutation_PI3K <- spread(mutation_PI3K,Cluster,MutatedSamples)
#Non  Low High
#I    II   III 
#99  136   56

mutation_PI3K$CIMP_H_percent <- mutation_PI3K$`CIMP-H`/56*100
mutation_PI3K$CIMP_L_percent <- mutation_PI3K$`CIMP-L`/136*100
mutation_PI3K$non_CIMP_percent <- mutation_PI3K$`non-CIMP`/99*100
write.table(mutation_PI3K,"mutation_PI3K.csv",sep=",")


#TGF-beta
TGF_beta_gene <- c("TGFBR1","TGFBR2","ACVR2A","ACVR1B","SMAD2",
                   "SMAD3","SMAD4","MYC")

mutation_I_TGF_beta <- mutation_I %>%
  filter(Hugo_Symbol %in% TGF_beta_gene )
mutation_III_TGF_beta <- mutation_III %>%
  filter(Hugo_Symbol %in% TGF_beta_gene )

mutation_II_TGF_beta <- mutation_II %>%
  filter(Hugo_Symbol %in% TGF_beta_gene )

#获得TGF_beta通路的突变数据
mutation_I_TGF_beta <- mutation_I_TGF_beta[,c(1,12)]
mutation_III_TGF_beta <- mutation_III_TGF_beta[,c(1,12)]
mutation_II_TGF_beta <- mutation_II_TGF_beta[,c(1,12)]


mutation_III_TGF_beta$Cluster <- "CIMP-H"
mutation_II_TGF_beta$Cluster <- "CIMP-L"
mutation_I_TGF_beta$Cluster <- "non-CIMP"
mutation_TGF_beta <- rbind(mutation_II_TGF_beta,mutation_III_TGF_beta,mutation_I_TGF_beta)
mutation_TGF_beta <- spread(mutation_TGF_beta,Cluster,MutatedSamples)

#Non  Low High
#I    II   III 
#99  136   56

mutation_TGF_beta$CIMP_H_percent <- mutation_TGF_beta$`CIMP-H`/56*100
mutation_TGF_beta$CIMP_L_percent <- mutation_TGF_beta$`CIMP-L`/136*100
mutation_TGF_beta$non_CIMP_percent <- mutation_TGF_beta$`non-CIMP`/99*100

write.table(mutation_TGF_beta,"mutation_TGF_beta.csv",sep=",")

#save("mutation_I","mutation_III","mutation_II",file = "mutation_29.RData")

#
##WNT通路
WNT_gene <- c("DKK1","DKK2","DKK3","DKK4","LRP5","FZD10","FAM123B",
              "AXIN2","APC","CTNNB1",
              "TCF7L2","FBXW7","ARID1A","SOX9")

#sum(mutation_III_WNT$Hugo_Symbol == "FAM123B")
mutation_I_WNT <- mutation_I %>%
  filter(Hugo_Symbol %in% WNT_gene )
mutation_III_WNT <- mutation_III %>%
  filter(Hugo_Symbol %in% WNT_gene )

mutation_II_WNT <- mutation_II %>%
  filter(Hugo_Symbol %in% WNT_gene)

#获得WNT通路的突变数据
mutation_I_WNT <- mutation_I_WNT[,c(1,12)]
mutation_III_WNT <- mutation_III_WNT[,c(1,12)]
mutation_II_WNT <- mutation_II_WNT[,c(1,12)]
########
mutation_III_WNT$Cluster <- "CIMP-H"
mutation_II_WNT$Cluster <- "CIMP-L"
mutation_I_WNT$Cluster <- "non-CIMP"

mutation_WNT <- rbind(mutation_II_WNT,mutation_III_WNT,mutation_I_WNT)
mutation_WNT <- spread(mutation_WNT,Cluster,MutatedSamples)


#Non  Low High
#I    II   III 
#99  136   56

mutation_WNT$CIMP_H_percent <- mutation_WNT$`CIMP-H`/56*100
mutation_WNT$CIMP_L_percent <- mutation_WNT$`CIMP-L`/136*100
mutation_WNT$non_CIMP_percent <- mutation_WNT$`non-CIMP`/99*100

write.table(mutation_WNT,"mutation_WNT.csv",sep=",")


#TP53通路
TP53_gene <- c("ATM","TP53")
mutation_I_TP53 <- mutation_I %>%
  filter(Hugo_Symbol %in% TP53_gene )
mutation_III_TP53 <- mutation_III %>%
  filter(Hugo_Symbol %in% TP53_gene )

mutation_II_TP53 <- mutation_II %>%
  filter(Hugo_Symbol %in% TP53_gene )

mutation_I_TP53 <- mutation_I_TP53[,c(1,12)]
mutation_III_TP53 <- mutation_III_TP53[,c(1,12)]
mutation_II_TP53 <- mutation_II_TP53[,c(1,12)]

mutation_III_TP53$Cluster <- "CIMP-H"
mutation_II_TP53$Cluster <- "CIMP-L"
mutation_I_TP53$Cluster <- "non-CIMP"


mutation_TP53<- rbind(mutation_II_TP53,mutation_III_TP53,mutation_I_TP53)
mutation_TP53 <- spread(mutation_TP53,Cluster,MutatedSamples)

#Non  Low High
#I    II   III 
#99  136   56

mutation_TP53$CIMP_H_percent <- mutation_TP53$`CIMP-H`/56*100
mutation_TP53$CIMP_L_percent <- mutation_TP53$`CIMP-L`/136*100
mutation_TP53$non_CIMP_percent <- mutation_TP53$`non-CIMP`/99*100

write.table(mutation_TP53,"mutation_TP53.csv",sep=",")
q()

