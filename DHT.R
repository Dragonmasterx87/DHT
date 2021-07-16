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
            features = c("GCG"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            #max.cutoff = 10,
            slot = 'data',
            order = TRUE)

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

# Some clusters including tuft, lymphocyte populations are missing
plot <- DimPlot(pancreas.integrated, reduction = "umap")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "schwann")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "lymphocytes")

# Saving this information in the metadata slot
table(Idents(pancreas.integrated))
pancreas.integrated$celltype <- Idents(pancreas.integrated)
head(pancreas.integrated@meta.data)

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "alpha", "delta", "gamma", 
               "ductal", "acinar", "stellate", 
               "schwann", "endothelium", "macrophage", "lymphocytes")
head(pancreas.integrated@meta.data$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated@meta.data$celltype <- factor(x = pancreas.integrated@meta.data$celltype, levels = my_levels)
Idents(pancreas.integrated) <- "celltype"
DimPlot(pancreas.integrated, split.by = "treatment", group.by = "celltype", label = TRUE)
DimPlot(pancreas.integrated, group.by = "treatment")
DimPlot(pancreas.integrated, split.by = "treatment", group.by = "DF.classifications_0.25_0.09_346",)

# Save file this will change but for showing them on 07132021 its fine
# saveRDS(pancreas.integrated, file = "D:/R-Projects/DHT/pancreas.integrated.rds")
pancreas.integrated <- readRDS("D:/R-Projects/DHT/pancreas.integrated.rds")

# Identify conserved cell markers
DefaultAssay(pancreas.integrated) <- "RNA"
markers.beta <- FindConservedMarkers(pancreas.integrated, ident.1 = "beta", grouping.var = "treatment", verbose = TRUE)
head(markers.beta)

markers.alpha <- FindConservedMarkers(pancreas.integrated, ident.1 = "alpha", grouping.var = "treatment", verbose = TRUE)
head(markers.alpha)

markers.delta <- FindConservedMarkers(pancreas.integrated, ident.1 = "delta", grouping.var = "treatment", verbose = TRUE)
head(markers.delta)

markers.gamma <- FindConservedMarkers(pancreas.integrated, ident.1 = "gamma", grouping.var = "treatment", verbose = TRUE)
head(markers.gamma)

markers.ductal <- FindConservedMarkers(pancreas.integrated, ident.1 = "ductal", grouping.var = "treatment", verbose = TRUE)
head(markers.ductal)

markers.acinar <- FindConservedMarkers(pancreas.integrated, ident.1 = "acinar", grouping.var = "treatment", verbose = TRUE)
head(markers.acinar)

markers.stellate <- FindConservedMarkers(pancreas.integrated, ident.1 = "stellate", grouping.var = "treatment", verbose = TRUE)
head(markers.stellate)

markers.schwann <- FindConservedMarkers(pancreas.integrated, ident.1 = "schwann", grouping.var = "treatment", verbose = TRUE)
head(markers.schwann)

markers.endothelium <- FindConservedMarkers(pancreas.integrated, ident.1 = "endothelium", grouping.var = "treatment", verbose = TRUE)
head(markers.endothelium)

markers.macrophage <- FindConservedMarkers(pancreas.integrated, ident.1 = "macrophage", grouping.var = "treatment", verbose = TRUE)
head(markers.macrophage)

markers.lymphocytes <- FindConservedMarkers(pancreas.integrated, ident.1 = "lymphocytes", grouping.var = "treatment", verbose = TRUE)
head(markers.lymphocytes)

# Identify conserved cell markers
Idents(pancreas.integrated) <- factor(Idents(pancreas.integrated), levels = c("beta", "alpha", "delta", "gamma", 
                                                                      "ductal", "acinar", "stellate", 
                                                                      "schwann", "endothelium", "macrophage", "lymphocytes"))
markers.to.plot <- c("INS", "MAFA", "IAPP", "GCG", "DPP4", "GC", "LEPR", "SST", "FRZB", "PPY", "CALB1", "THSD7A",
                     "CFTR", "TFPI2", "MMP7", "CELA2A", "CELA2B", "CELA3A", "COL3A1", "FMOD", "PDGFRB", 
                     "SOX10", "CDH19", "NGFR", "CD34", "ENG", "VWF", "CD86", "CSF1R", "FCER1G", "NKG7", "IL2RB", "CCL5")

# Advanced coding for ggplot2
# Create a new metadata slot containing combined info, segregating clusters and samples
Idents(object = pancreas.integrated) <- "celltype"
pancreas.integrated$celltype.sample <- paste(Idents(pancreas.integrated),pancreas.integrated$treatment, sep = "_")
table(pancreas.integrated@meta.data$celltype.sample)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("beta_EtOH", "beta_DHT[10nM]", "alpha_EtOH", "alpha_DHT[10nM]", "delta_EtOH", "delta_DHT[10nM]", 
                "gamma_EtOH", "gamma_DHT[10nM]", "ductal_EtOH", "ductal_DHT[10nM]", "acinar_EtOH", "acinar_DHT[10nM]", 
                "stellate_EtOH", "stellate_DHT[10nM]", "schwann_EtOH", "schwann_DHT[10nM]", "endothelium_EtOH", "endothelium_DHT[10nM]", 
                "macrophage_EtOH", "macrophage_DHT[10nM]", "lymphocytes_EtOH", "lymphocytes_DHT[10nM]")
head(pancreas.integrated@meta.data$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated@meta.data$celltype.sample <- factor(x = pancreas.integrated@meta.data$celltype.sample, levels = my_levels2)
table(pancreas.integrated@meta.data$celltype.sample)

# Re select organized idents
Idents(pancreas.integrated) <- "celltype.sample"
DotPlot(pancreas.integrated,  
        dot.scale = 8, 
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient(low =c("white"), high =c("red")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

# Older dotplot configuration, shows only percentage not expression
Idents(pancreas.integrated) <- "celltype"
DotPlot(pancreas.integrated, features = rev(markers.to.plot), 
        cols = c("blue", "red"), 
        dot.scale = 8, 
        split.by = "treatment") + 
  RotatedAxis() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  theme_light() + 
  #coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =10, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =8, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10))

# Diff gene testing across conditions
pancreas.integrated$treatment.dht <- paste(Idents(pancreas.integrated), pancreas.integrated$treatment, sep = "_")
pancreas.integrated$celltype.split <- Idents(pancreas.integrated)
Idents(pancreas.integrated) <- "treatment.dht"
beta.DHT.response <- FindMarkers(pancreas.integrated, ident.1 = "beta_EtOH", ident.2 = "beta_DHT[10nM]", verbose = FALSE)
head(beta.DHT.response, n = 15)


plots <- VlnPlot(beta.cells, features = c("MT-CO3", "MAF", "AR"), group.by = "treatment", 
                 pt.size = 0, combine = TRUE)
wrap_plots(plots = plots, nrow = 1)



VlnPlot(pancreas.integrated, features = c("MT-CO3", "MT-ND1", "MT-ATP6", "CA2", "PDK4"), group.by = "treatment")
write.csv(beta.DHT.response, file = "D:/R-Projects/DHT/Data output/beta.DHT.response.csv")


X1 <- NULL
table(x = FetchData(pancreas.integrated, vars = c('celltype', 'sex')))
x1 <- subset(pancreas.integrated, subset = (celltype == c("beta", "alpha")) & (sex == "Male"))
table(x = FetchData(x1, vars = c('celltype', 'sex')))
x2 <- subset(pancreas.integrated, subset = (celltype != c("beta")) & (sex != "Male")) # wont run because you cant subset a vector with no cells which is what is left once all male cells are removed :)
table(x = FetchData(x2, vars = c('celltype', 'sex')))






theme_set(theme_cowplot())
beta.cells <- subset(pancreas.integrated, idents = "beta")
Idents(beta.cells) <- "treatment"
avg.beta.cells <- log1p(AverageExpression(beta.cells, verbose = FALSE)$RNA)
avg.beta.cells$gene <- rownames(avg.beta.cells)

alpha.cells <- subset(pancreas.integrated, idents = "alpha")
Idents(alpha.cells) <- "treatment"
avg.alpha.cells <- log1p(AverageExpression(alpha.cells, verbose = FALSE)$RNA)
avg.alpha.cells$gene <- rownames(avg.alpha.cells)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(avg.beta.cells, aes("ctrl", "DHT[10nM]")) + geom_point() + ggtitle("Beta Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
plot_grid(p1, p2)



# Extra
# Installation of the latest released version
install.packages('GOplot')
library(GOplot)
packageVersion("GOplot")

gene.data <- read.csv("D:/R-Projects/DHT/Data output/beta.DHT.de.data.csv")
up.go <- read.csv("D:/R-Projects/DHT/Data output/up/DHTGOup.csv")
down.go <- read.csv("D:/R-Projects/DHT/Data output/Down/DHTGOdown.csv")

head(gene.data)
head(up.go)
circ <- circle_dat(up.go, gene.data)
circ
process <- List('cellular response to decreased oxygen levels', "cellular response to hypoxia",
                'response to unfolded protein', 'canonical glycolysis',
                'glucose catabolic process to pyruvate', 'glycolytic process through glucose-6-phosphate',
                'gluconeogenesis', 'cellular response to oxidative stress',
                'protein stabilization ', 'amino acid transport ')

process
chord <- chord_dat(data = circ, genes = gene.data)
chord <- chord_dat(data = circ, process = process)
chord <- chord_dat(data = circ, genes = gene.data, process = process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)


GOBubble(circ, labels = 1)


