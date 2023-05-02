library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(harmony)
library(SeuratData)
library(SeuratDisk)

# Download Haniffa's skin scRNA-seq data from https:/app.cellatlas.io/diseased-skin/dataset/3/download
Convert("submission.h5ad", dest = "h5seurat", overwrite = TRUE)
healthy <- LoadH5Seurat("./submission.h5seurat")
writeRDS(healthy,"seurat.RDS")

healthy<-readRDS("seurat.RDS")

healthy$batch<-healthy$sample_id
DefaultAssay(healthy)<-"RNA"
healthy <- FindVariableFeatures(healthy, selection.method = "vst", nfeatures = 3000)

dir.create("./HS_combined_v2/")
setwd("./HS_combined_v2/")

HS22_Myeloid_S4.data <- Read10X(data.dir = "/share/studies/Dermatology_Data/HS Data Portal/GSE155850/HS22_Myeloid_S4/outs/filtered_feature_bc_matrix/")
HS22_Myeloid_S4b.data <- Read10X(data.dir = "/share/studies/Dermatology_Data/HS Data Portal/GSE155850/HS22_Myeloid_S4b/outs/filtered_feature_bc_matrix/")

intersect_HS22_Myeloid_S4<-intersect(unlist(HS22_Myeloid_S4.data@Dimnames[2]),unlist(HS22_Myeloid_S4b.data@Dimnames[2]))
HS22_Myeloid_S4.data1<-HS22_Myeloid_S4.data[,intersect_HS22_Myeloid_S4]
HS22_Myeloid_S4b.data1<-HS22_Myeloid_S4b.data[,intersect_HS22_Myeloid_S4]
HS22_Myeloid_S4all.data<-HS22_Myeloid_S4.data1+HS22_Myeloid_S4b.data1

HS22_Myeloid_S4 <- CreateSeuratObject(counts = HS22_Myeloid_S4all.data, project = "HS22_Myeloid_S4", min.cells = 0, min.features = 0)
HS22_Myeloid_S4 <- RenameCells(object = HS22_Myeloid_S4, add.cell.id = "HS22_Myeloid_S4")

HS22_Myeloid_S4[["percent.mt"]] <- PercentageFeatureSet(HS22_Myeloid_S4, pattern = "^MT-")
pdf("HS22_Myeloid_S4_QC_plot_before_trimming.pdf", width = 16, height = 5)
VlnPlot(HS22_Myeloid_S4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
HS22_Myeloid_S4 <- subset(HS22_Myeloid_S4, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 20)
pdf("HS22_Myeloid_S4_QC_plot_after_trimming.pdf", width = 16, height = 5)
VlnPlot(HS22_Myeloid_S4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

HSEx17Myeloid_S8.data <- Read10X(data.dir = "/share/studies/Dermatology_Data/HS Data Portal/GSE155850/HSEx17Myeloid_S8/outs/filtered_feature_bc_matrix/")
HSEx17Myeloid_S8b.data <- Read10X(data.dir = "/share/studies/Dermatology_Data/HS Data Portal/GSE155850/HSEx17Myeloid_S8b/outs/filtered_feature_bc_matrix/")

intersect_HSEx17Myeloid_S8<-intersect(unlist(HSEx17Myeloid_S8.data@Dimnames[2]),unlist(HSEx17Myeloid_S8b.data@Dimnames[2]))
HSEx17Myeloid_S8.data1<-HSEx17Myeloid_S8.data[,intersect_HSEx17Myeloid_S8]
HSEx17Myeloid_S8b.data1<-HSEx17Myeloid_S8b.data[,intersect_HSEx17Myeloid_S8]
HSEx17Myeloid_S8all.data<-HSEx17Myeloid_S8.data1+HSEx17Myeloid_S8b.data1

HSEx17Myeloid_S8 <- CreateSeuratObject(counts = HSEx17Myeloid_S8all.data, project = "HSEx17Myeloid_S8", min.cells = 0, min.features = 0)
HSEx17Myeloid_S8 <- RenameCells(object = HSEx17Myeloid_S8, add.cell.id = "HSEx17Myeloid_S8")

HSEx17Myeloid_S8[["percent.mt"]] <- PercentageFeatureSet(HSEx17Myeloid_S8, pattern = "^MT-")
pdf("HSEx17Myeloid_S8_QC_plot_before_trimming.pdf", width = 16, height = 5)
VlnPlot(HSEx17Myeloid_S8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
HSEx17Myeloid_S8 <- subset(HSEx17Myeloid_S8, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 20)
pdf("HSEx17Myeloid_S8_QC_plot_after_trimming.pdf", width = 16, height = 5)
VlnPlot(HSEx17Myeloid_S8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

HSDerm190115_S4.data <- Read10X(data.dir = "/share/studies/Dermatology_Data/HS Data Portal/GSE154775/HSDerm190115_S4/outs/filtered_feature_bc_matrix/")
HSDerm190115_S4 <- CreateSeuratObject(counts = HSDerm190115_S4.data, project = "HSDerm190115_S4", min.cells = 0, min.features = 0)
HSDerm190115_S4 <- RenameCells(object = HSDerm190115_S4, add.cell.id = "HSDerm190115_S4")

HSDerm190115_S4[["percent.mt"]] <- PercentageFeatureSet(HSDerm190115_S4, pattern = "^MT-")
pdf("HSDerm190115_S4_QC_plot_before_trimming.pdf", width = 16, height = 5)
VlnPlot(HSDerm190115_S4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
HSDerm190115_S4 <- subset(HSDerm190115_S4, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 20)
pdf("HSDerm190115_S4_QC_plot_after_trimming.pdf", width = 16, height = 5)
VlnPlot(HSDerm190115_S4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

HSEpi190115_S3.data <- Read10X(data.dir = "/share/studies/Dermatology_Data/HS Data Portal/GSE154775/HSEpi190115_S3/outs/filtered_feature_bc_matrix/")
HSEpi190115_S3 <- CreateSeuratObject(counts = HSEpi190115_S3.data, project = "HSEpi190115_S3", min.cells = 0, min.features = 0)
HSEpi190115_S3 <- RenameCells(object = HSEpi190115_S3, add.cell.id = "HSEpi190115_S3")

HSEpi190115_S3[["percent.mt"]] <- PercentageFeatureSet(HSEpi190115_S3, pattern = "^MT-")
pdf("HSEpi190115_S3_QC_plot_before_trimming.pdf", width = 16, height = 5)
VlnPlot(HSEpi190115_S3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
HSEpi190115_S3 <- subset(HSEpi190115_S3, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 20)
pdf("HSEpi190115_S3_QC_plot_after_trimming.pdf", width = 16, height = 5)
VlnPlot(HSEpi190115_S3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

HSDerm190103_S2.data <- Read10X(data.dir = "/share/studies/Dermatology_Data/HS Data Portal/GSE154775/HSDerm190103_S2/outs/filtered_feature_bc_matrix/")
HSDerm190103_S2 <- CreateSeuratObject(counts = HSDerm190103_S2.data, project = "HSDerm190103_S2", min.cells = 0, min.features = 0)
HSDerm190103_S2 <- RenameCells(object = HSDerm190103_S2, add.cell.id = "HSDerm190103_S2")

HSDerm190103_S2[["percent.mt"]] <- PercentageFeatureSet(HSDerm190103_S2, pattern = "^MT-")
pdf("HSDerm190103_S2_QC_plot_before_trimming.pdf", width = 16, height = 5)
VlnPlot(HSDerm190103_S2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
HSDerm190103_S2 <- subset(HSDerm190103_S2, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 20)
pdf("HSDerm190103_S2_QC_plot_after_trimming.pdf", width = 16, height = 5)
VlnPlot(HSDerm190103_S2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

HSEpi190103_S1.data <- Read10X(data.dir = "/share/studies/Dermatology_Data/HS Data Portal/GSE154775/HSEpi190103_S1/outs/filtered_feature_bc_matrix/")
HSEpi190103_S1 <- CreateSeuratObject(counts = HSEpi190103_S1.data, project = "HSEpi190103_S1", min.cells = 0, min.features = 0)
HSEpi190103_S1 <- RenameCells(object = HSEpi190103_S1, add.cell.id = "HSEpi190103_S1")

HSEpi190103_S1[["percent.mt"]] <- PercentageFeatureSet(HSEpi190103_S1, pattern = "^MT-")
pdf("HSEpi190103_S1_QC_plot_before_trimming.pdf", width = 16, height = 5)
VlnPlot(HSEpi190103_S1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
HSEpi190103_S1 <- subset(HSEpi190103_S1, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 20)
pdf("HSEpi190103_S1_QC_plot_after_trimming.pdf", width = 16, height = 5)
VlnPlot(HSEpi190103_S1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

HS_lesional_dermis_S2.data <- Read10X(data.dir = "/share/studies/Dermatology_Data/HS Data Portal/GSE154775/HS_lesional_dermis_S2/outs/filtered_feature_bc_matrix/")
HS_lesional_dermis_S2 <- CreateSeuratObject(counts = HS_lesional_dermis_S2.data, project = "HS_lesional_dermis_S2", min.cells = 0, min.features = 0)
HS_lesional_dermis_S2 <- RenameCells(object = HS_lesional_dermis_S2, add.cell.id = "HS_lesional_dermis_S2")

HS_lesional_dermis_S2[["percent.mt"]] <- PercentageFeatureSet(HS_lesional_dermis_S2, pattern = "^MT-")
pdf("HS_lesional_dermis_S2_QC_plot_before_trimming.pdf", width = 16, height = 5)
VlnPlot(HS_lesional_dermis_S2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
HS_lesional_dermis_S2 <- subset(HS_lesional_dermis_S2, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 20)
pdf("HS_lesional_dermis_S2_QC_plot_after_trimming.pdf", width = 16, height = 5)
VlnPlot(HS_lesional_dermis_S2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

HS_lesional_epi_S1.data <- Read10X(data.dir = "/share/studies/Dermatology_Data/HS Data Portal/GSE154775/HS_lesional_epi_S1/outs/filtered_feature_bc_matrix/")
HS_lesional_epi_S1 <- CreateSeuratObject(counts = HS_lesional_epi_S1.data, project = "HS_lesional_epi_S1", min.cells = 0, min.features = 0)
HS_lesional_epi_S1 <- RenameCells(object = HS_lesional_epi_S1, add.cell.id = "HS_lesional_epi_S1")

HS_lesional_epi_S1[["percent.mt"]] <- PercentageFeatureSet(HS_lesional_epi_S1, pattern = "^MT-")
pdf("HS_lesional_epi_S1_QC_plot_before_trimming.pdf", width = 16, height = 5)
VlnPlot(HS_lesional_epi_S1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
HS_lesional_epi_S1 <- subset(HS_lesional_epi_S1, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 20)
pdf("HS_lesional_epi_S1_QC_plot_after_trimming.pdf", width = 16, height = 5)
VlnPlot(HS_lesional_epi_S1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

HS2_S3.data <- Read10X(data.dir = "/share/studies/Dermatology_Data/HS Data Portal/GSE154775/HS2_S3/outs/filtered_feature_bc_matrix/")
HS2_S3 <- CreateSeuratObject(counts = HS2_S3.data, project = "HS2_S3", min.cells = 0, min.features = 0)
HS2_S3 <- RenameCells(object = HS2_S3, add.cell.id = "HS2_S3")

HS2_S3[["percent.mt"]] <- PercentageFeatureSet(HS2_S3, pattern = "^MT-")
pdf("HS2_S3_QC_plot_before_trimming.pdf", width = 16, height = 5)
VlnPlot(HS2_S3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
HS2_S3 <- subset(HS2_S3, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 20)
pdf("HS2_S3_QC_plot_after_trimming.pdf", width = 16, height = 5)
VlnPlot(HS2_S3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

HS1_S4.data <- Read10X(data.dir = "/share/studies/Dermatology_Data/HS Data Portal/GSE154775/HS1_S4/outs/filtered_feature_bc_matrix/")
HS1_S4 <- CreateSeuratObject(counts = HS1_S4.data, project = "HS1_S4", min.cells = 0, min.features = 0)
HS1_S4 <- RenameCells(object = HS1_S4, add.cell.id = "HS1_S4")

HS1_S4[["percent.mt"]] <- PercentageFeatureSet(HS1_S4, pattern = "^MT-")
pdf("HS1_S4_QC_plot_before_trimming.pdf", width = 16, height = 5)
VlnPlot(HS1_S4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
HS1_S4 <- subset(HS1_S4, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 20)
pdf("HS1_S4_QC_plot_after_trimming.pdf", width = 16, height = 5)
VlnPlot(HS1_S4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

HS_07.data <- Read10X(data.dir = "/share/studies/Dermatology_Data/Haniffa_Dataset/data_download/HS_07/outs/filtered_feature_bc_matrix/")
HS_07 <- CreateSeuratObject(counts = HS_07.data, project = "HS_07", min.cells = 0, min.features = 0)
HS_07 <- RenameCells(object = HS_07, add.cell.id = "HS_07")

HS_07[["percent.mt"]] <- PercentageFeatureSet(HS_07, pattern = "^MT-")
pdf("HS_07_QC_plot_before_trimming.pdf", width = 16, height = 5)
VlnPlot(HS_07, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
HS_07 <- subset(HS_07, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 20)
pdf("HS_07_QC_plot_after_trimming.pdf", width = 16, height = 5)
VlnPlot(HS_07, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

hs_subset<-merge(HS22_Myeloid_S4,y=c(HSEx17Myeloid_S8,HSDerm190115_S4,HSEpi190115_S3,HSDerm190103_S2,HSEpi190103_S1,HS_lesional_dermis_S2,HS_lesional_epi_S1,HS2_S3,HS1_S4,HS_07))

hs_subset$batch<-hs_subset$orig.ident
saveRDS(hs_subset,"hs.rds")
DefaultAssay(hs_subset)<-"RNA"
hs_subset <- NormalizeData(hs_subset)
hs_subset <- FindVariableFeatures(hs_subset, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(hs_subset)
hs_subset <- ScaleData(hs_subset, features = all.genes, vars.to.regress = "percent.mt")
hs_subset <- RunPCA(hs_subset, features = VariableFeatures(object = hs_subset))
hs_subset <- RunHarmony(hs_subset, "batch", plot_convergence = TRUE, assay.use ="RNA")

hs_subset <- FindNeighbors(hs_subset, dims = 1:30, verbose = TRUE, reduction = "harmony")
hs_subset <- FindClusters(hs_subset, verbose = TRUE, resolution = 1)
hs_subset <- RunUMAP(hs_subset, dims = 1:30, verbose = TRUE, reduction = "harmony")
saveRDS(hs_subset, "hs_subset.rds")
hs_subset<-readRDS("hs_subset.rds")
transfer.anchors <- FindTransferAnchors(reference = healthy, query = hs_subset, features = VariableFeatures(object = healthy),
    reference.assay = "RNA", query.assay = "RNA", reduction = "pcaproject")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = healthy$full_clustering,
    weight.reduction = hs_subset[["harmony"]], dims = 1:30)

hs_subset <- AddMetaData(hs_subset, metadata = celltype.predictions)
saveRDS(hs_subset, "hs_subset.rds")

hs_subset$predicted.id2<-hs_subset$predicted.id
hs_subset$predicted.id2[hs_subset@meta.data$prediction.score.max<=0.5]<-"Unknown"

saveRDS(hs_subset, "hs_subset.rds")

dir.create("./resolution_1.5/")
setwd("./resolution_1.5/")

hs_subset<-readRDS("../hs_subset.rds")

hs_subset <- FindClusters(hs_subset, verbose = TRUE, resolution = 1.5)
saveRDS(hs_subset, "hs_subset.rds")
hs_subset<-readRDS("hs_subset.rds")

dir.create("./HS_only_resolution_1.5_annot/")
setwd("./HS_only_resolution_1.5_annot/")

hs_subset$clusters_orig<-hs_subset@active.ident
hs_subset$clusters_annot1<-as.character(hs_subset@active.ident)

hs_subset$clusters_annot1[hs_subset$clusters_orig==0 ]<-"Basal KCs"
hs_subset$clusters_annot1[hs_subset$clusters_orig==1 ]<-"ILCs/NK Cells"
hs_subset$clusters_annot1[hs_subset$clusters_orig==2 ]<-"T Cells"
hs_subset$clusters_annot1[hs_subset$clusters_orig==3 ]<-"Macrophages"
hs_subset$clusters_annot1[hs_subset$clusters_orig==4 ]<-"Basal-Stem KCs"
hs_subset$clusters_annot1[hs_subset$clusters_orig==5 ]<-"Basal KCs"
hs_subset$clusters_annot1[hs_subset$clusters_orig==6 ]<-"Suprabasal KCs"
hs_subset$clusters_annot1[hs_subset$clusters_orig==7 ]<-"C1q-hi Macrophages"
hs_subset$clusters_annot1[hs_subset$clusters_orig==8 ]<-"B Cells"
hs_subset$clusters_annot1[hs_subset$clusters_orig==9 ]<-"Basal-Stem KCs"
hs_subset$clusters_annot1[hs_subset$clusters_orig==10 ]<-"Monocytes"
hs_subset$clusters_annot1[hs_subset$clusters_orig==11 ]<-"T Cells"
hs_subset$clusters_annot1[hs_subset$clusters_orig==12 ]<-"Suprabasal KCs"
hs_subset$clusters_annot1[hs_subset$clusters_orig==13 ]<-"Endothelial Cells"
hs_subset$clusters_annot1[hs_subset$clusters_orig==14 ]<-"Proliferative KCs"
hs_subset$clusters_annot1[hs_subset$clusters_orig==15 ]<-"Pericytes"
hs_subset$clusters_annot1[hs_subset$clusters_orig==16 ]<-"Fibroblasts"
hs_subset$clusters_annot1[hs_subset$clusters_orig==17 ]<-"Langerhans Cells"
hs_subset$clusters_annot1[hs_subset$clusters_orig==18 ]<-"Suprabasal KCs"
hs_subset$clusters_annot1[hs_subset$clusters_orig==19 ]<-"Mast Cells"
hs_subset$clusters_annot1[hs_subset$clusters_orig==20 ]<-"Plasma"
hs_subset$clusters_annot1[hs_subset$clusters_orig==21 ]<-"Proliferative KCs"
hs_subset$clusters_annot1[hs_subset$clusters_orig==22 ]<-"Melanocytes"
hs_subset$clusters_annot1[hs_subset$clusters_orig==23 ]<-"Proliferative KCs"
hs_subset$clusters_annot1[hs_subset$clusters_orig==24 ]<-"Basal KCs"
hs_subset$clusters_annot1[hs_subset$clusters_orig==25 ]<-"Basal-Stem KCs"
hs_subset$clusters_annot1[hs_subset$clusters_orig==26 ]<-"Basal KCs"
hs_subset$clusters_annot1[hs_subset$clusters_orig==27 ]<-"Resident Memory T Cells"
hs_subset$clusters_annot1[hs_subset$clusters_orig==28 ]<-"Dendritic Cells"
hs_subset$clusters_annot1[hs_subset$clusters_orig==29 ]<-"ILCs/NK Cells"
hs_subset$clusters_annot1[hs_subset$clusters_orig==30 ]<-"ILCs/NK Cells"
hs_subset$clusters_annot1<-as.factor(hs_subset$clusters_annot1)
hs_subset@active.ident<-hs_subset$clusters_annot1

saveRDS(hs_subset, "hs_subset.rds")
hs_subset.slim <- DietSeurat(hs_subset, counts = TRUE, data = TRUE, scale.data = FALSE, dimreducs = c('umap'))
saveRDS(hs_subset.slim, "HS_merged_lite.RDS")

pdf("DimPlot_UMAP_by_ManualAnnot_LabelT.pdf", width = 10, height = 6)
DimPlot(hs_subset, reduction = "umap", group.by="clusters_annot1",label=T,raster=F,cols=cols)
dev.off()
pdf("DimPlot_UMAP_by_ManualAnnot_LabelF.pdf", width = 10, height = 6)
DimPlot(hs_subset, reduction = "umap", group.by="clusters_annot1",label=F,raster=F,cols=cols)
dev.off()
pdf("DimPlot_UMAP_by_FindClusters.pdf", width = 8, height = 6)
DimPlot(hs_subset, reduction = "umap",label=T,raster=F,cols=cols,group.by="clusters_orig")
dev.off()
pdf("DimPlot_UMAP_by_Samples.pdf", width = 8, height = 6)
DimPlot(hs_subset, reduction = "umap",label=T,raster=F,group.by="batch",cols=cols)
dev.off()


DefaultAssay(hs_subset)<-"RNA"

pdf("Freq_barplot.pdf", width = 12, height = 6)
xx=barplot(table(hs_subset@active.ident))
xx
dat_freq<-as.numeric(table(hs_subset@active.ident))
text(x = xx, y = dat_freq, label =dat_freq, pos = 1, cex = 0.8)
dev.off()
dat_freq

hs_subset.markers <- FindAllMarkers(hs_subset, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(hs_subset.markers,"hs_subset.markers.markers_all.csv",sep=",",row.names=F)

hs_subset.markers <- FindAllMarkers(hs_subset, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(hs_subset.markers,"hs.markers.markers_positive.csv",sep=",", row.names=F)

#hs_subset.markers<-read.table("hs_subset.markers.markers_positive.csv", sep=",",header=T)

top10 <-hs_subset.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
write.table(top10, "hs.markers.top10_genes.csv", sep = ',', row.names=F)

top2 <- hs_subset.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)

top30 <- hs_subset.markers %>% group_by(cluster) %>% top_n(30, avg_log2FC)
write.table(top30, "hs.markers.top30_genes.csv", sep = ',', row.names=F)

pdf("DoHeatmap_genes.pdf", width = 24, height = 20)
DoHeatmap(object = hs_subset, features = as.character(top10$gene))
dev.off()

library(ggplot2)
plot1 <- DotPlot(object = hs_subset, features = unique(as.character(top10$gene)), dot.scale=5,scale=T)
plot1_new<- plot1 + theme(axis.text.x = element_text(angle = 90))

pdf("BubblePlot_genes.pdf", width = 50, height = 10)
plot1_new
dev.off()

DefaultAssay(hs_subset)<-"RNA"

gene_list<-as.character(top30$gene)
unlink("./UMAP_plots_top30_genes/", recursive = TRUE)
dir.create("./UMAP_plots_top30_genes/")
for (i in 1:length(gene_list)){
  file_png_name<-paste("./UMAP_plots_top30_genes/UMAP_",gene_list[i],".pdf",sep="")
  fplot<-FeaturePlot(object = hs_subset, features = gene_list[i],reduction = 'umap') + scale_color_viridis()
  pdf(file_png_name, width=4.5, height=4)
  print(fplot)
  dev.off()
}
unlink("./Violin_plots_top30_genes/", recursive = TRUE)
dir.create("./Violin_plots_top30_genes/")
for (i in 1:length(gene_list)){
  file_png_name<-paste("./Violin_plots_top30_genes/VlnPlot_",gene_list[i],".pdf",sep="")
  fplot<-VlnPlot(object = hs_subset, features = gene_list[i], pt.size = 0, combine = FALSE)
  pdf(file_png_name, width=4.5, height=4)
  print(fplot)
  dev.off()
}

DefaultAssay(hs_subset)<-"RNA"
gene_list<-c("CD19","CD3D","CD3E","CD3G","FCGR3A","FCGR3B","ITGAX","NCAM1","CD14","CD8A","PTPRC","KLRB1","CD38","CD4","KRT8","KRT18","ITGA6","KRT15","TP63","PPP3CA","CAV2","CTNNAL1","COL17A1","CDK1","PCNA","MTX1","KRT1","STMN1","ICAM1","CD83","CCL20","TNF","BDKRB1","KRT19","CNFN","KRT10")


dir.create("./UMAP_plots_select_genes/")
for (i in 1:length(gene_list)){
  if (paste0(gene_list[i]) %in% hs_subset@assays$RNA@counts@Dimnames[[1]] == 'FALSE')next()
  file_png_name<-paste("./UMAP_plots_select_genes/UMAP_",gene_list[i],".pdf",sep="")
  fplot<-FeaturePlot(object = hs_subset, features = gene_list[i], reduction = 'umap') + scale_color_viridis()

  pdf(file_png_name, width=4.5, height=4)

  print(fplot)
  dev.off()
}
dir.create("./Violin_plots_select_genes/")
for (i in 1:length(gene_list)){
  if (paste0(gene_list[i]) %in% hs_subset@assays$RNA@counts@Dimnames[[1]] == 'FALSE')next()
  file_png_name<-paste("./Violin_plots_select_genes/VlnPlot_",gene_list[i],".pdf",sep="")
  fplot<-VlnPlot(object = hs_subset, combine = FALSE, features = gene_list[i], pt.size = 0)
  pdf(file_png_name, width=4.5, height=4)
  print(fplot)
  dev.off()
}

