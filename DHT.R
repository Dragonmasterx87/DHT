# LOAD LIBRARIES ####
# Restart Rstudio or R
# Run the following code once you have Seurat installed
suppressWarnings({
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
  library(DoubletFinder)})

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle3")

# Set global environment parameter
options(future.globals.maxSize = 8000 * 1024^2)

# OBJECT SETUP AND NORMALIZATION ####
# STEP 1: Load 10X data ####
HP2107001_ctrl.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Raw Data\F7a GEX_HP21070_01\filtered_feature_bc_matrix)")
HP2107001_DHT.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Raw Data\F7b GEX_HP21070_01_DHT\filtered_feature_bc_matrix)")
HP2107701_ctrl.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Raw Data\F8a_GEX_HP21077_01\filtered_feature_bc_matrix)")
HP2107701_DHT.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Raw Data\F8b_GEX_HP21077_01_DHT\filtered_feature_bc_matrix)")
HP2107901_ctrl.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Raw Data\F9a_GEX_HP21079_01\filtered_feature_bc_matrix)")
HP2107901_DHT.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Raw Data\F9b_GEX_HP21079_01_DHT\filtered_feature_bc_matrix)")

# STEP 2: Create Seurat objects ####
HP2107001_ctrl <- CreateSeuratObject(counts = HP2107001_ctrl.data, min.features = 500)
HP2107001_DHT <- CreateSeuratObject(counts = HP2107001_DHT.data, min.features = 500)
HP2107701_ctrl <- CreateSeuratObject(counts = HP2107701_ctrl.data, min.features = 500)
HP2107701_DHT <- CreateSeuratObject(counts = HP2107701_DHT.data, min.features = 500)
HP2107901_ctrl <- CreateSeuratObject(counts = HP2107901_ctrl.data, min.features = 500)
HP2107901_DHT <- CreateSeuratObject(counts = HP2107901_DHT.data, min.features = 500)

# Sample specific Metadata addition
HP2107001_ctrl$sample <- "HP2107001_ctrl"
HP2107001_DHT$sample <- "HP2107001_DHT"
HP2107701_ctrl$sample <- "HP2107701_ctrl"
HP2107701_DHT$sample <- "HP2107701_DHT"
HP2107901_ctrl$sample <- "HP2107901_ctrl"
HP2107901_DHT$sample <- "HP2107901_DHT"

# Sex specific Metadata addition
HP2107001_ctrl$sex <- "Male"
HP2107001_DHT$sex <- "Male"
HP2107701_ctrl$sex <- "Male"
HP2107701_DHT$sex <- "Male"
HP2107901_ctrl$sex <- "Male"
HP2107901_DHT$sex <- "Male"

# Treatment specific Metadata addition
HP2107001_ctrl$treatment <- "EtOH"
HP2107001_DHT$treatment <- "DHT[10nM]"
HP2107701_ctrl$treatment <- "EtOH"
HP2107701_DHT$treatment <- "DHT[10nM]"
HP2107901_ctrl$treatment <- "EtOH"
HP2107901_DHT$treatment <- "DHT[10nM]"

# STEP 3: Thresholding ####
# The operator can add columns to object metadata. This is a great place to stash QC stats
HP2107001_ctrl[["percent.mt"]] <- PercentageFeatureSet(object = HP2107001_ctrl, pattern = "^MT-")
HP2107001_DHT[["percent.mt"]] <- PercentageFeatureSet(object = HP2107001_DHT, pattern = "^MT-")
HP2107701_ctrl[["percent.mt"]] <- PercentageFeatureSet(object = HP2107701_ctrl, pattern = "^MT-")
HP2107701_DHT[["percent.mt"]] <- PercentageFeatureSet(object = HP2107701_DHT, pattern = "^MT-")
HP2107901_ctrl[["percent.mt"]] <- PercentageFeatureSet(object = HP2107901_ctrl, pattern = "^MT-")
HP2107901_DHT[["percent.mt"]] <- PercentageFeatureSet(object = HP2107901_DHT, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(object = HP2107001_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107001_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107701_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107701_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107901_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107901_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# RNA based cell thresholding
HP2107001_ctrl <- subset(x = HP2107001_ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HP2107001_DHT <- subset(x = HP2107001_DHT, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HP2107701_ctrl <- subset(x = HP2107701_ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HP2107701_DHT <- subset(x = HP2107701_DHT, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HP2107901_ctrl <- subset(x = HP2107901_ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
HP2107901_DHT <- subset(x = HP2107901_DHT, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)

# Visualize QC metrics post thresholding as a violin plot
VlnPlot(object = HP2107001_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107001_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107701_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107701_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107901_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107901_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Step 4: Add cell IDs ####
# Add cell IDs
HP2107001_ctrl <- RenameCells(HP2107001_ctrl, add.cell.id = "HP2107001_ctrl")
HP2107001_DHT <- RenameCells(HP2107001_DHT, add.cell.id = "HP2107001_DHT")
HP2107701_ctrl <- RenameCells(HP2107701_ctrl, add.cell.id = "HP2107701_ctrl")
HP2107701_DHT <- RenameCells(HP2107701_DHT, add.cell.id = "HP2107701_DHT")
HP2107901_ctrl <- RenameCells(HP2107901_ctrl, add.cell.id = "HP2107901_ctrl")
HP2107901_DHT <- RenameCells(HP2107901_DHT, add.cell.id = "HP2107901_DHT")

# Step 5: Merge Datasets
# Based on comment to Issue #4753 https://github.com/satijalab/seurat/issues/4753
# We use RPCA to yield conserved mapping and set Tx as control refrence samples
# Merge panc_sex datasets
pancreas.list <- list("HP2107001_ctrl" = HP2107001_ctrl, "HP2107001_DHT" = HP2107001_DHT,
                      "HP2107701_ctrl" = HP2107701_ctrl, "HP2107701_DHT" = HP2107701_DHT,
                      "HP2107901_ctrl" = HP2107901_ctrl, "HP2107901_DHT" = HP2107901_DHT)

# Step 6: Data normalization
#Normalise data
pancreas.list <- lapply(X = pancreas.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

# Step 7: Feature selection
# Select features for downstream integration
pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list)
pancreas.list <- lapply(X = pancreas.list, FUN = function(x) {
  x <- ScaleData(x, features = pancreas.features, verbose = FALSE)
  x <- RunPCA(x, features = pancreas.features, verbose = FALSE)
})

# Step 8: Anchor identification and data integration
# Identify anchors and integrate dataset
pancreas.list[c(1, 3, 5)] # check that you are correctly picking up control datasets for refrence integration
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, reference = c(1,3,5), # Takes 18min 13sec to run when using cca as reduction
                                           reduction = "rpca", dims = 1:30, verbose = TRUE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30, verbose = TRUE)

# Step 9: Linear dimensionality assessment
# Look at your default assay
DefaultAssay(object = pancreas.integrated)

# Change default assay to integrated, to view dimensionality
DefaultAssay(object = pancreas.integrated) <- "integrated"

# Scaling this is weird, but as done in https://satijalab.org/seurat/articles/integration_large_datasets.html
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)

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

# Based of clustree assessment choose res = 0.3, conservative approach
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
Idents(pancreas.integrated) <- "sex"
Idents(pancreas.integrated) <- "sample"
Idents(pancreas.integrated) <- "seurat_clusters"
DimPlot(pancreas.integrated, reduction = "umap", label = TRUE)

#Visualize gene expression
DefaultAssay(object = pancreas.integrated) <- "RNA"
DefaultAssay(object = pancreas.integrated)
FeaturePlot(object = pancreas.integrated,
            features = c("MKI67"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            #min.cutoff = 0,
            #max.cutoff = 10,
            slot = 'counts',
            order = TRUE)

#Rename Idents
pancreas.integrated <- RenameIdents(pancreas.integrated, 
                                    "0" = "Alpha", 
                                    "1" = "Beta",
                                    "2" = "Ductal", 
                                    "3" = "Beta",
                                    "4" = "Transdifferentiating Endocrine", 
                                    "5" = "Acinar",
                                    "6" = "Activated Stellate", 
                                    "7" = "Quiescent Stellate",
                                    "8" = "Endothelial", 
                                    "9" = "Delta",
                                    "10" = "Alpha", 
                                    "11" = "Ductal",
                                    "12" = "Gamma",
                                    "13" = "Ductal",
                                    "14" = "Macrophage",
                                    "15" = "Proliferating Stellate",
                                    "16" = "Schwann",
                                    "17" = "Mast",
                                    "18" = "T-Lymphocyte"
                                    )

DimPlot(pancreas.integrated, reduction = "umap", label = TRUE)

# Saving this information in the metadata slot
table(Idents(pancreas.integrated))
pancreas.integrated$celltype <- Idents(pancreas.integrated)
head(pancreas.integrated@meta.data)

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("Beta", "Alpha", "Transdifferentiating Endocrine", "Delta", "Gamma", 
               "Ductal", "Acinar", 
               "Quiescent Stellate", "Activated Stellate", "Proliferating Stellate",
               "Macrophage", "T-Lymphocyte", "Mast",
               "Schwann", "Endothelial")
head(pancreas.integrated@meta.data$celltype)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated@meta.data$celltype <- factor(x = pancreas.integrated@meta.data$celltype, levels = my_levels)
Idents(pancreas.integrated) <- "celltype"

# Observing cells
DimPlot(pancreas.integrated, split.by = "sample", group.by = "celltype", label = FALSE, ncol = 2)
DimPlot(pancreas.integrated, group.by = "treatment")
DimPlot(pancreas.integrated, reduction = "umap", 
        cols = c("darkgreen",
                 "red",
                 "sienna3",
                 "mediumseagreen",
                 "turquoise4",
                 "black",
                 "royalblue1",
                 "yellow4",
                 "gray30",
                 "darkred",
                 "orange2",
                 "darkmagenta",
                 "deeppink2"
                 ),
                 label = TRUE)

# Save file this will change but for showing them on 07132021 its fine
saveRDS(pancreas.integrated, file = "C:/Users/mqadir/Box/Lab 2301/RNAseq DHT data/wkdir/pancreas.integrated.rds")
pancreas.integrated <- readRDS("C:/Users/mqadir/Box/Lab 2301/RNAseq DHT data/wkdir/pancreas.integrated.rds")

# Identify conserved cell markers
DefaultAssay(pancreas.integrated) <- "RNA"
markers.beta <- FindConservedMarkers(pancreas.integrated, ident.1 = "Beta", grouping.var = "treatment", verbose = TRUE)
head(markers.beta)

markers.alpha <- FindConservedMarkers(pancreas.integrated, ident.1 = "Alpha", grouping.var = "treatment", verbose = TRUE)
head(markers.alpha)

markers.tranendo <- FindConservedMarkers(pancreas.integrated, ident.1 = "Transdifferentiating Endocrine", grouping.var = "treatment", verbose = TRUE)
head(markers.tranendo)

markers.delta <- FindConservedMarkers(pancreas.integrated, ident.1 = "Delta", grouping.var = "treatment", verbose = TRUE)
head(markers.delta)

markers.gamma <- FindConservedMarkers(pancreas.integrated, ident.1 = "Gamma", grouping.var = "treatment", verbose = TRUE)
head(markers.gamma)

markers.ductal <- FindConservedMarkers(pancreas.integrated, ident.1 = "Ductal", grouping.var = "treatment", verbose = TRUE)
head(markers.ductal)

markers.acinar <- FindConservedMarkers(pancreas.integrated, ident.1 = "Acinar", grouping.var = "treatment", verbose = TRUE)
head(markers.acinar)

markers.quiescentstellate <- FindConservedMarkers(pancreas.integrated, ident.1 = "Quiescent Stellate", grouping.var = "treatment", verbose = TRUE)
head(markers.quiescentstellate)

markers.activatedstellate <- FindConservedMarkers(pancreas.integrated, ident.1 = "Activated Stellate", grouping.var = "treatment", verbose = TRUE)
head(markers.activatedstellate)

markers.prolifstellate <- FindConservedMarkers(pancreas.integrated, ident.1 = "Proliferating Stellate", grouping.var = "treatment", verbose = TRUE)
head(markers.prolifstellate)

markers.macrophage <- FindConservedMarkers(pancreas.integrated, ident.1 = "Macrophage", grouping.var = "treatment", verbose = TRUE)
head(markers.macrophage)

markers.tlympho <- FindConservedMarkers(pancreas.integrated, ident.1 = "T-Lymphocyte", grouping.var = "treatment", verbose = TRUE)
head(markers.tlympho)

markers.mast <- FindConservedMarkers(pancreas.integrated, ident.1 = "Mast", grouping.var = "treatment", verbose = TRUE)
head(markers.mast)

markers.schwann <- FindConservedMarkers(pancreas.integrated, ident.1 = "Schwann", grouping.var = "treatment", verbose = TRUE)
head(markers.schwann)

markers.endothelial <- FindConservedMarkers(pancreas.integrated, ident.1 = "Endothelial", grouping.var = "treatment", verbose = TRUE)
head(markers.endothelial)

# Identify conserved cell markers
Idents(pancreas.integrated) <- factor(Idents(pancreas.integrated), levels = c("Beta", "Alpha", "Transdifferentiating Endocrine", "Delta", "Gamma", 
                                                                              "Ductal", "Acinar", 
                                                                              "Quiescent Stellate", "Activated Stellate", "Proliferating Stellate", 
                                                                              "Macrophage", "T-Lymphocyte", "Mast", "Schwann", "Endothelial"))
markers.to.plot <- c("INS", "IAPP", "NKX6-1", "MAFA", "MAFB", "GCG", "DPP4", "GC", "LEPR", "SST", "FRZB", "PPY", "CALB1", "THSD7A",
                     "CFTR", "KRT19", "MMP7", "CELA2A", "CELA2B", "CELA3A", "COL3A1", "FMOD", "PDGFRB", "MKI67", "HIST1H4C", "STMN1", 
                     "CD86", "CSF1R", "SDS", "NKG7", "IL2RB", "CCL5", "RGS13", "TPSB2", "TPSAB1", "SOX10", "CDH19", "NGFR", "CD34", "ENG", "VWF")

# Advanced coding for ggplot2
# Create a new metadata slot containing combined info, segregating clusters and samples
Idents(object = pancreas.integrated) <- "celltype"
pancreas.integrated$celltype.sample <- paste(Idents(pancreas.integrated),pancreas.integrated$treatment, sep = "_")
table(pancreas.integrated@meta.data$celltype.sample)

# New metadata column is not paired, so we need to pair
my_levels2 <- c("Beta_EtOH", "Beta_DHT[10nM]", "Transdifferentiating Endocrine_EtOH", "Transdifferentiating Endocrine_DHT[10nM]", "Alpha_EtOH", "Alpha_DHT[10nM]", "Delta_EtOH", "Delta_DHT[10nM]", "Gamma_EtOH", "Gamma_DHT[10nM]", 
                "Ductal_EtOH", "Ductal_DHT[10nM]", "Acinar_EtOH", "Acinar_DHT[10nM]", 
                "Quiescent Stellate_EtOH", "Quiescent Stellate_DHT[10nM]", "Activated Stellate_EtOH", "Activated Stellate_DHT[10nM]", "Proliferating Stellate_EtOH", "Proliferating Stellate_DHT[10nM]",
                "Macrophage_EtOH", "Macrophage_DHT[10nM]", "T-Lymphocyte_EtOH", "T-Lymphocyte_DHT[10nM]", "Mast_EtOH", "Mast_DHT[10nM]", "Schwann_EtOH", "Schwann_DHT[10nM]", "Endothelial_EtOH", "Endothelial_DHT[10nM]")
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
# pancreas.integrated$treatment.dht <- paste(Idents(pancreas.integrated), pancreas.integrated$treatment, sep = "_")
# pancreas.integrated$celltype.split <- Idents(pancreas.integrated)
Idents(pancreas.integrated) <- "celltype.sample"
beta.DHT.response <- FindMarkers(pancreas.integrated, 
                                 ident.1 = "Beta_DHT[10nM]", ident.2 = "Beta_EtOH", 
                                 test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                 min.pct = 0.1,
                                 logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                 pseudocount.use = 1,
                                 verbose = FALSE)
head(beta.DHT.response, n = 15)
write.csv(beta.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\beta.DHT.response.csv)")

alpha.DHT.response <- FindMarkers(pancreas.integrated, 
                                 ident.1 = "Alpha_DHT[10nM]", ident.2 = "Alpha_EtOH", 
                                 test.use = "wilcox", # Based on #2938 DESeq2 not recommended for single cell gene expression analysis
                                 min.pct = 0.1,
                                 logfc.threshold = 0.137504, # based on output log2 so 0.137504 is ~1.1 FC
                                 pseudocount.use = 1,
                                 verbose = FALSE)
head(alpha.DHT.response, n = 15)
write.csv(alpha.DHT.response, file = r"(C:\Users\mqadir\Box\Lab 2301\RNAseq DHT data\Data output\alpha.DHT.response.csv)")

# Plotting DE genes
Idents(pancreas.integrated) <- "celltype"
beta.cells <- subset(pancreas.integrated, idents = "Beta")
plots <- VlnPlot(beta.cells, features = c("INS", "DDIT3", "MIF", "DEPP1", "PLCG2", "IAPP"), group.by = "treatment", 
                 pt.size = 1, combine = TRUE)
plots <- VlnPlot(beta.cells, features = c("MT-CO3", "MT-ND1", "MT-ND4", "MT-ATP6", "MT-CO1", "MT-CYB"), group.by = "treatment", 
                 pt.size = 1, combine = TRUE)
wrap_plots(plots = plots, nrow = 1, ncol = 1)


# Calculating percentages
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


