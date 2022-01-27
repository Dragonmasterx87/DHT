# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 
# R version 4.0.2 (2020-06-22) -- "Taking Off Again"

# SESSION INFORMATION AND SETUP ####
sessionInfo()
R.Version()

# Set working directory or save as a project in a specific folder in your PC
setwd(r"(C:\Users\mqadir\Box\!FAHD\10. GLP1-E2 Project\R_wkdir)")

# Check for Working Directory
getwd()


# pACKAGE INSTALLATION ####
# Package instalation
# Install biocmanager
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

# Install computational packages
# Biocmanager dependant
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
BiocManager::install("Glimma")
BiocManager::install("apeglm")

# Directly from CRAN
install.packages(c("R.basic"), contriburl="http://www.braju.com/R/repos/")
install.packages("survival")
install.packages("digest")
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("calibrate")
install.packages("gplots")

# load Packages
suppressWarnings(
  {
    library(limma)
    library(edgeR)
    library(Glimma)
    library(gplots)
    library(RColorBrewer)
    library(Matrix)
    library(survival)
    library(digest)
    library(tidyr)
    library(tidyverse)
    library(pheatmap)
    library(calibrate)
    library(pheatmap)
    }
  )

# LOADING DATA ####
# Read data into R
# For RPPA data. Data is read as a DF you need to change to matrix for heatmap.2
MS.heatmap <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Metabolomics\metabdata_AA.csv)", header = TRUE, sep = ",", row.names = 1)
MS.heatmap <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Metabolomics\metabdata_GTO.csv)", header = TRUE, sep = ",", row.names = 1)
MS.heatmap <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Metabolomics\metabdata_PC.csv)", header = TRUE, sep = ",", row.names = 1)
MS.heatmap <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Metabolomics\metabdata_F4.csv)", header = TRUE, sep = ",", row.names = 1)

MS.heatmap

# Convert to matrix
MS.heatmap.matrix <- data.matrix(MS.heatmap)

# DATA ANALYSIS
# As defined by Kevin https://www.biostars.org/p/284721/#284723 and https://www.biostars.org/p/318628/
# Based on https://www.biostars.org/p/318628/
# Color scheme
require("RColorBrewer")
myCol <- colorRampPalette(c("dodgerblue4", "white", "red3"))(100)
myBreaks <- seq(-2, 2, length.out=101)

#Transform to Z-scale
heat <- t(scale(t(MS.heatmap)))

# Plot
require("gplots")

#Euclidean distance
heatmap.2(heat,
          #Rowv=as.dendrogram(hr),
          #Colv=as.dendrogram(hc),
          col=myCol,
          breaks=myBreaks,
          main="Phosphate Carrier",
          key=T, keysize=1.0,
          scale="none",
          density.info="none",
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
          trace="none",
          sepwidth = c(0.001,0.001),
          sepcolor = "black",
          colsep=1:ncol(heat),
          rowsep=1:nrow(heat),
          cexRow=1.0, #size of gene font
          cexCol=1.2, #size of sample font
          distfun=function(x) dist(x, method="euclidean"),
          hclustfun=function(x) hclust(x, method="ward.D2"),
          margins = c(10, 12))

# create heatmap using pheatmap
pheatmap(heat,
         color = colorRampPalette(c("dodgerblue4", "white", "red3"))(100),
         cluster_cols = FALSE, #clustering columns
         cluster_rows = FALSE, #clustering rows
         gaps_col = c(1, 4, 7), #gaps between columns
         border_color = "black" #color of borders between genes/proteins
         )

# DATA ANALYSIS: HEATMAPS ####
# Hierarchical clustering with heatmaps
## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlGn")
morecols <- colorRampPalette(mypalette)

# Set up colour vector for celltype variable
# my_palette <- colorRampPalette(c("darkgreen", "yellow", "red"))(n = 50)
# col.cell <- c('lightblue', 'blue', 'red')

mypalette <- brewer.pal(11,"RdYlGn")
morecols <- colorRampPalette(mypalette)     


heatmap.2(MS.heatmap.matrix,
          col= rev(morecols(50)),
          trace="none", 
          main="RPPA",
          #ColSideColors=col.cell,
          #dendrogram = 'row',
          Rowv = TRUE,
          Colv = TRUE,
          #key = NULL,
          key.title = "Z-Score",
          #ColSideColors=col.cell,
          scale="row",
          sepwidth=c(0.1,0.1),
          #rowsep = c(2, 5, 9, 10, 12, 15, 17, 20, 21),
          margins = c(8, 11),
          breaks = seq(-10,10, length.out = 50),
          cexRow = c(1))


