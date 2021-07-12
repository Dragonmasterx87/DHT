# LOAD LIBRARIES ####
# Restart Rstudio or R
# Run the following code once you have Seurat installed
library(ggplot2)
library(cowplot)
library(Matrix)
library(ggridges)
library(ggrepel)
library(dplyr)
library(Seurat)
library(monocle3)
library(plotly)
library(clustree)
library(patchwork)
library(future)
library(DoubletFinder)

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle3")

# Set global environment parameter
options(future.globals.maxSize = 4000 * 1024^2)

# OBJECT SETUP AND NORMALIZATION ####
# STEP 1: Load 10X data ####
HP2107001_ctrl.data <- Read10X(data.dir = "D:/R-Projects/DHT/F7a GEX_HP21070_01/filtered_feature_bc_matrix/")
HP2107001_DHT.data <- Read10X(data.dir = "D:/R-Projects/DHT/F7b GEX_HP21070_01_DHT/filtered_feature_bc_matrix/")

# STEP 2: Create Seurat objects ####
HP2107001_ctrl <- CreateSeuratObject(counts = HP2107001_ctrl.data)
HP2107001_DHT <- CreateSeuratObject(counts = HP2107001_DHT.data)

# Sample specific Metadata addition
HP2107001_ctrl$sample <- "HP2107001_ctrl"
HP2107001_ctrl$sex <- "Male"
HP2107001_ctrl$treatment <- "EtOH"
HP2107001_DHT$sample <- "HP2107001_DHT"
HP2107001_DHT$sex <- "Male"
HP2107001_DHT$treatment <- "DHT[10nM]"

# STEP 3: Thresholding ####
# The operator can add columns to object metadata. This is a great place to stash QC stats
HP2107001_ctrl[["percent.mt"]] <- PercentageFeatureSet(object = HP2107001_ctrl, pattern = "^MT-")
HP2107001_DHT[["percent.mt"]] <- PercentageFeatureSet(object = HP2107001_DHT, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(object = HP2107001_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107001_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# RNA based cell thresholding
HP2107001_ctrl <- subset(x = HP2107001_ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HP2107001_DHT <- subset(x = HP2107001_DHT, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)

# Visualize QC metrics post thresholding as a violin plot
VlnPlot(object = HP2107001_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107001_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Step 4: Add cell IDs ####
# Add cell IDs
HP2107001_ctrl <- RenameCells(HP2107001_ctrl, add.cell.id = "HP2107001_ctrl")
HP2107001_DHT <- RenameCells(HP2107001_DHT, add.cell.id = "HP2107001_DHT")

# Step 5: Merge Datasets
# Merge panc_sex datasets
pancreas.list <- list("HP2107001_ctrl" = HP2107001_ctrl, "HP2107001_DHT" = HP2107001_DHT)

# Step 6: Data normalization
# Normalize the dataset using SCTransform
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], verbose = TRUE)
}

# Step 7: Feature selection
# Select features for downstream integration
pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features,
                                    verbose = TRUE)

# Step 8: Anchor identification and data integration
# Identify anchors and integrate dataset
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT",
                                           anchor.features = pancreas.features, verbose = TRUE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT",
                                     verbose = TRUE)

# Step 9: Linear dimensionality assessment
# Look at your default assay
DefaultAssay(object = pancreas.integrated)

# Change default assay to integrated, to view dimensionality
DefaultAssay(object = pancreas.integrated) <- "integrated"

# Dimensionality assessment using PCA analysis
pancreas.integrated <- RunPCA(pancreas.integrated, features = VariableFeatures(object = pancreas.integrated))

# Examine data dimensionality
ElbowPlot(pancreas.integrated)

# Step 9a: CLUSTER ANALYSIS ####
# Clustering, we run multiple permutations to allow clustree to analyze optimal clustering resolution.
pancreas.integrated <- FindNeighbors(object = pancreas.integrated, dims = 1:30)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.1)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.2)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.3)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.4)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.5)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.6)

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution
clustree(pancreas.integrated, prefix = "integrated_snn_res.")

# Based of clustree assessment choose res = 0.3, stable yet biologically relevant
# Beyond 0.4 massive cluster destabilization occurs
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.3)

# Alternatively build a cluster tree
DefaultAssay(object = pancreas.integrated) <- "integrated"
pancreas.integrated=BuildClusterTree(pancreas.integrated, slot = "scale.data")
PlotClusterTree(pancreas.integrated)

# Step 10: non-linear dimensionality assessment ####
# Run PCA and UMAP calculations
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30)

# Change default assay to integrated, to view dimensionality
Idents(pancreas.integrated) <- "treatment"
DimPlot(pancreas.integrated, reduction = "umap", label = FALSE)

#Visualize gene expression
DefaultAssay(object = pancreas.integrated) <- "RNA"
FeaturePlot(object = pancreas.integrated,
            features = c("VWF"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            #max.cutoff = 10,
            slot = 'data',
            order = TRUE)

# At this point, I discovered that there are a cluster of cells which contain a mixture of alpha-beta cells
# This is possibly doublets which are contaminating the analysis.

# Discoverng doublets for ctrl
ctrl <- HP2107001_ctrl

## Pre-process Seurat object (standard)
DefaultAssay(object = ctrl) <- "RNA"
ctrl <- NormalizeData(ctrl)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)
ctrl <- ScaleData(ctrl)
ctrl <- RunPCA(ctrl)
ctrl <- RunUMAP(ctrl, dims = 1:10)
DimPlot(ctrl, reduction = "umap", label = FALSE)

## Clustering
ctrl <- FindNeighbors(object = ctrl, dims = 1:10)
ctrl <- FindClusters(object = ctrl, resolution = 0.3)
DimPlot(ctrl, reduction = "umap", label = FALSE)

## pK Identification (no ground-truth)
sweep.res.list_ctrl <- paramSweep_v3(ctrl, PCs = 1:10, sct = FALSE)
sweep.stats_ctrl <- summarizeSweep(sweep.res.list_ctrl, GT = FALSE)
bcmvn_ctrl <- find.pK(sweep.stats_ctrl)

## Homotypic Doublet Proportion Estimate
annotations <- ctrl$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  ## ex: annotations <- ctrl@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(ctrl@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
ctrl <- doubletFinder_v3(ctrl, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(ctrl, reduction = "umap", label = FALSE, group.by = "DF.classifications_0.25_0.09_346")
# singlets_HP2107001_ctrl <- subset(ctrl, subset = DF.classifications_0.25_0.09_353 == "Singlet")
# singlets_HP2107001_ctrl
# HP2107001_ctrl

## Discoverng doublets for DHT
DHT <- HP2107001_DHT

## Pre-process Seurat object (standard)
DefaultAssay(object = DHT) <- "RNA"
DHT <- NormalizeData(DHT)
DHT <- FindVariableFeatures(DHT, selection.method = "vst", nfeatures = 2000)
DHT <- ScaleData(DHT)
DHT <- RunPCA(DHT)
DHT <- RunUMAP(DHT, dims = 1:10)
DimPlot(DHT, reduction = "umap", label = FALSE)

## Clustering
DHT <- FindNeighbors(object = DHT, dims = 1:10)
DHT <- FindClusters(object = DHT, resolution = 0.3)
DimPlot(DHT, reduction = "umap", label = FALSE)

## pK Identification (no ground-truth)
sweep.res.list_DHT <- paramSweep_v3(DHT, PCs = 1:10, sct = FALSE)
sweep.stats_DHT <- summarizeSweep(sweep.res.list_DHT, GT = FALSE)
bcmvn_DHT <- find.pK(sweep.stats_DHT)

## Homotypic Doublet Proportion Estimate
annotations <- DHT$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  ## ex: annotations <- DHT@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(DHT@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
DHT <- doubletFinder_v3(DHT, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(DHT, reduction = "umap", label = FALSE, group.by = "DF.classifications_0.25_0.09_488")
# singlets_HP2107001_DHT <- subset(DHT, subset = DF.classifications_0.25_0.09_494 == "Singlet")
# singlets_HP2107001_DHT
# HP2107001_DHT

## Extract metadata
# ctrl_singlet_meta <- ctrl$DF.classifications_0.25_0.09_353
# DHT_singlet_meta <- DHT$DF.classifications_0.25_0.09_494
# doublets = c(ctrl_singlet_meta, DHT_singlet_meta) # Concatenate https://www.tutorialspoint.com/how-to-concatenate-two-or-more-vectors-in-r 
# there is also another way using pylr but I havent looked into it to much yet: https://github.com/satijalab/seurat/issues/2081
# pancreas.integrated$doublet <- doublets 
# table(pancreas.integrated$doublet)

# Subsetting out doublets
# pancreas.integrated.singlets <- subset(pancreas.integrated, subset = doublet == "Singlet")
# pancreas.integrated
# pancreas.integrated.singlets
# DimPlot(pancreas.integrated, reduction = "umap", group.by = "doublet", label = FALSE)
# DimPlot(pancreas.integrated.singlets, reduction = "umap", group.by = "doublet", label = FALSE)

# One will have to subset out singlets from individual data first
HP2107001_ctrl_singlets <- subset(ctrl, subset = DF.classifications_0.25_0.09_346 == "Singlet")
HP2107001_DHT_singlets <- subset(DHT, subset = DF.classifications_0.25_0.09_488 == "Singlet")
HP2107001_ctrl_singlets
HP2107001_ctrl
HP2107001_DHT_singlets
HP2107001_DHT

# Lets try this again
# OBJECT SETUP AND NORMALIZATION ####
# STEP 1: Load 10X data ####
# STEP 2: Create Seurat objects ####
HP2107001_ctrl <- HP2107001_ctrl_singlets
HP2107001_DHT <- HP2107001_DHT_singlets

# Sample specific Metadata addition
HP2107001_ctrl$sample <- "HP2107001_ctrl"
HP2107001_ctrl$sex <- "Male"
HP2107001_ctrl$treatment <- "EtOH"
HP2107001_DHT$sample <- "HP2107001_DHT"
HP2107001_DHT$sex <- "Male"
HP2107001_DHT$treatment <- "DHT[10nM]"

# STEP 3: Thresholding ####
# The operator can add columns to object metadata. This is a great place to stash QC stats
HP2107001_ctrl[["percent.mt"]] <- PercentageFeatureSet(object = HP2107001_ctrl, pattern = "^MT-")
HP2107001_DHT[["percent.mt"]] <- PercentageFeatureSet(object = HP2107001_DHT, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(object = HP2107001_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107001_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# RNA based cell thresholding
HP2107001_ctrl <- subset(x = HP2107001_ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HP2107001_DHT <- subset(x = HP2107001_DHT, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)

# Visualize QC metrics post thresholding as a violin plot
VlnPlot(object = HP2107001_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107001_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Step 4: Add cell IDs ####
# Add cell IDs
HP2107001_ctrl <- RenameCells(HP2107001_ctrl, add.cell.id = "HP2107001_ctrl")
HP2107001_DHT <- RenameCells(HP2107001_DHT, add.cell.id = "HP2107001_DHT")

# Step 5: Merge Datasets
# Merge panc_sex datasets
pancreas.list <- list("HP2107001_ctrl" = HP2107001_ctrl, "HP2107001_DHT" = HP2107001_DHT)

# Step 6: Data normalization
# Normalize the dataset using SCTransform
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], verbose = TRUE)
}

# Step 7: Feature selection
# Select features for downstream integration
pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features,
                                    verbose = TRUE)

# Step 8: Anchor identification and data integration
# Identify anchors and integrate dataset
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT",
                                           anchor.features = pancreas.features, verbose = TRUE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT",
                                     verbose = TRUE)

# Step 9: Linear dimensionality assessment
# Look at your default assay
DefaultAssay(object = pancreas.integrated)

# Change default assay to integrated, to view dimensionality
DefaultAssay(object = pancreas.integrated) <- "integrated"

# Dimensionality assessment using PCA analysis
pancreas.integrated <- RunPCA(pancreas.integrated, features = VariableFeatures(object = pancreas.integrated))

# Examine data dimensionality
ElbowPlot(pancreas.integrated)

# Step 9a: CLUSTER ANALYSIS ####
# Clustering, we run multiple permutations to allow clustree to analyze optimal clustering resolution.
pancreas.integrated <- FindNeighbors(object = pancreas.integrated, dims = 1:20)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.1)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.2)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.3)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.4)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.5)
pancreas.integrated <- FindClusters(object = pancreas.integrated, resolution = 0.6)

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution
clustree(pancreas.integrated, prefix = "integrated_snn_res.")

# Based of clustree assessment choose res = 0.3, stable yet biologically relevant
# Beyond 0.4 massive cluster destabilization occurs
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.5)

# Alternatively build a cluster tree
DefaultAssay(object = pancreas.integrated) <- "integrated"
pancreas.integrated=BuildClusterTree(pancreas.integrated, slot = "scale.data")
PlotClusterTree(pancreas.integrated)

# Step 10: non-linear dimensionality assessment ####
# Run PCA and UMAP calculations
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:40) #this will change with the remaining n

# Change default assay to integrated, to view dimensionality
Idents(pancreas.integrated) <- "treatment"
DimPlot(pancreas.integrated, reduction = "umap", label = TRUE, group.by = "seurat_clusters")

#Visualize gene expression
DefaultAssay(object = pancreas.integrated) <- "RNA"
FeaturePlot(object = pancreas.integrated,
            features = c("INS", "GCG"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            #max.cutoff = 10,
            slot = 'data',
            order = TRUE)

# Save file this will change but for showing them on 07132021 its fine
# saveRDS(pancreas.integrated, file = "D:/R-Projects/DHT/pancreas.integrated.rds")

#New metadata column
pancreas.integrated$new_clusters <- pancreas.integrated$integrated_snn_res.0.5
Idents(pancreas.integrated) <- "new_clusters"
DimPlot(pancreas.integrated, reduction = "umap", label = TRUE)

#Rename Idents
pancreas.integrated <- RenameIdents(pancreas.integrated, 
                                    "0" = "beta", 
                                    "1" = "alpha",
                                    "2" = "ductal", 
                                    "3" = "acinar",
                                    "4" = "beta", 
                                    "5" = "beta",
                                    "6" = "alpha", 
                                    "7" = "endothelium",
                                    "8" = "stellate", 
                                    "9" = "beta",
                                    "10" = "beta", 
                                    "11" = "alpha",
                                    "12" = "gamma",
                                    "13" = "macrophage",
                                    "14" = "delta",
                                    "15" = "stellate"
                                    )
DimPlot(pancreas.integrated, reduction = "umap", label = TRUE)
