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
setwd(r"(C:\Users\mqadir\Box\Lab 2301\Coding Scripts\RPPA)")

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
    }
  )

# LOADING DATA ####
# Read data into R
MS.heatmap <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\Coding Scripts\Mass Spec\mass_spec_actual.csv)", header = TRUE, sep = ",", row.names = 1)
MS.heatmap <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\Coding Scripts\Mass Spec\mass_spec_actual_DHT.csv)", header = TRUE, sep = ",", row.names = 1)
MS.heatmap <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\Coding Scripts\Mass Spec\mass_spec_actual_DHT_GLP1.csv)", header = TRUE, sep = ",", row.names = 1)
MS.heatmap <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\Coding Scripts\Mass Spec\mass_spec_actual_Negatively_corelated.csv)", header = TRUE, sep = ",", row.names = 1)

# For RPPA
MS.heatmap <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\Coding Scripts\Mass Spec\Selected total proteins_for heatmap.csv)", header = TRUE, sep = ",", row.names = 1)

# DATA ANALYSIS: HEATMAPS ####
# Hierarchical clustering with heatmaps
## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlGn")
morecols <- colorRampPalette(mypalette)

# Set up colour vector for celltype variable
my_palette <- colorRampPalette(c("darkblue", "yellow", "darkred"))(n = 299)
col.cell <- c('lightblue', 'blue', 'red'
              )

MS.heatmap.matrix <- data.matrix(MS.heatmap)
heatmap.2(MS.heatmap.matrix,
          col=rev(morecols(50)),
          trace="none", 
          main="RPPA",
          #ColSideColors=col.cell,
          dendrogram = 'row',
          Rowv = FALSE,
          Colv = FALSE,
          #key = NULL,
          key.title = "Z-Score",
          scale="row",
          sepwidth=c(0.2,0.2),
          #rowsep = c(2, 5, 9, 10, 12, 15, 17, 20, 21),
          margins = c(3, 15),
          cexRow = c(1))

