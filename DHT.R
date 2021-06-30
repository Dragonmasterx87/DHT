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

#################### #
#################### #

# STEP 3: Thresholding ####
# The operator can add columns to object metadata. This is a great place to stash QC stats
HP2107001_ctrl[["percent.mt"]] <- PercentageFeatureSet(object = HP2107001_ctrl, pattern = "^MT-")
HP2107001_DHT[["percent.mt"]] <- PercentageFeatureSet(object = HP2107001_DHT, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(object = HP2107001_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107001_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# RNA based cell thresholding
HP2107001_ctrl <- subset(x = HP2107001_ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25)
HP2107001_DHT <- subset(x = HP2107001_DHT, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25)

# Visualize QC metrics post thresholding as a violin plot
VlnPlot(object = HP2107001_ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = HP2107001_DHT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#################### #
#################### #

# Step 4: Add cell IDs ####
# Add cell IDs
HP2107001_ctrl <- RenameCells(HP2107001_ctrl, add.cell.id = "HP2107001_ctrl")
HP2107001_DHT <- RenameCells(HP2107001_DHT, add.cell.id = "HP2107001_DHT")

#################### #
#################### #

# Step 5: Merge Datasets
# Merge panc_sex datasets
pancreas.list <- list("HP2107001_DHT" = HP2107001_DHT, "HP2107001_ctrl" = HP2107001_ctrl)

# Step 6: Data normalization
# Normalize the dataset using SCTransform
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], verbose = TRUE)
}

#################### #
#################### #

# Step 7: Feature selection
# Select features for downstream integration
pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features, 
                                    verbose = TRUE)
#################### #
#################### #

# Step 8: Anchor identification and data integration
# Identify anchors and integrate dataset
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = TRUE)

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", 
                                     verbose = TRUE)

#################### #
#################### #

# Step 9: Linear dimensionality assessment
# Look at your default assay
DefaultAssay(object = pancreas.integrated)

# Change default assay to integrated, to view dimensionality
DefaultAssay(object = pancreas.integrated) <- "integrated"

# Dimensionality assessment using PCA analysis
pancreas.integrated <- RunPCA(pancreas.integrated, features = VariableFeatures(object = pancreas.integrated))

# Examine data dimensionality
ElbowPlot(pancreas.integrated)

#################### #
#################### #

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
# We will return to cluster analysis after optimally UMAP cluster info

# Step 10: non-linear dimensionality assessment ####
# Run PCA and UMAP calculations
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30)

# Change default assay to integrated, to view dimensionality
Idents(pancreas.integrated) <- "treatment"
DimPlot(pancreas.integrated, reduction = "umap", label = FALSE)

#Visualize gene expression
DefaultAssay(object = pancreas.integrated) <- "RNA"
FeaturePlot(object = pancreas.integrated, 
            features = c("INS", "GCG", "SST", "PPY", "GHRL"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            #max.cutoff = 10,
            slot = 'data',
            order = TRUE)
