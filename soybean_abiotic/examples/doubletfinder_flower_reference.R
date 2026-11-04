# Source: 图1.docx, extracted text lines 1-39
# 原文草莓花双细胞示例；非大豆分析入口
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

##参考下面去除双细胞QC
#!/usr/bin/env R
library(Seurat)
library(dplyr)
library(DoubletFinder)

folders<-c("Flower_ck1","Flower_ck2","Flower_inf_20h1","Flower_inf_20h2","Flower_mock_20hpt1","Flower_mock_20hpt2")
for (i in folders){
  a<-CreateSeuratObject(counts = Read10X(paste0("/public/home/houzy/sc-seq/Fve_Flowers/Fve_flower0.2/Seuratdata/",i),gene.column = 1),project = i, min.cells = 3, min.features = 200)
  a<-RenameCells(a, add.cell.id =i)
  a@meta.data$Group<-c(rep("Fve",ncol(a)))
  a <-NormalizeData(a) %>% 
    ScaleData(verbose = FALSE) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    RunPCA(verbose = FALSE) %>% 
    RunTSNE(dims=1:30) %>% 
    RunUMAP(dims=1:30) %>%
    FindNeighbors( dims = 1:30) %>%
    FindClusters(resolution = 0.7)
  #### Doubletfinder
  sweep.res.list <- paramSweep(a, PCs = 1:30, sct = F,num.cores = 52)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  ## 排除不能检出的同源doublets，优化期望的doublets数量
  annotations <- a@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)  
  DoubletRate = ncol(a)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计???
  #估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
  nExp_poi <- round(DoubletRate*length(a$seurat_clusters))  #最好提供celltype，而不是seurat_clusters???
  # 计算双细胞比???
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  a <- doubletFinder(a, PCs = 1:30, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
  # 存储seurat对象
  assign(paste0('Fve_',i),a)
}


save(Fve_Flower_ck1,Fve_Flower_ck2,Fve_Flower_inf_20h1,Fve_Flower_inf_20h2,Fve_Flower_mock_20hpt1,Fve_Flower_mock_20hpt2,file = "/public/home/houzy/sc-seq/Fve_Flowers/Fve_flower0.2/Seurat/01.QC/Fve_Flower.Rdata")
