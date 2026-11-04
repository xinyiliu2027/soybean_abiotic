# Source: Figure 1 document, extracted text lines 1-39
# Strawberry flower doublet example from the source; not a soybean analysis entry point.
# See README.md and docs/KNOWN_LIMITATIONS.md before running.

# Reference workflow for doublet removal and quality control.
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
  # Adjust the expected doublet count for undetectable homotypic doublets.
  annotations <- a@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)  
  DoubletRate = ncol(a)*8*1e-6 # Assumes an increase of 0.8 percentage points per 1,000 cells; verify this assumption.
  # Estimate the homotypic doublet proportion using seurat_clusters in modelHomotypic().
  nExp_poi <- round(DoubletRate*length(a$seurat_clusters))  # Consider whether validated cell-type annotations should replace seurat_clusters.
  # Adjust the expected doublet count; verify the underlying rate assumption.
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  a <- doubletFinder(a, PCs = 1:30, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
  # Store the Seurat object.
  assign(paste0('Fve_',i),a)
}


save(Fve_Flower_ck1,Fve_Flower_ck2,Fve_Flower_inf_20h1,Fve_Flower_inf_20h2,Fve_Flower_mock_20hpt1,Fve_Flower_mock_20hpt2,file = "/public/home/houzy/sc-seq/Fve_Flowers/Fve_flower0.2/Seurat/01.QC/Fve_Flower.Rdata")
