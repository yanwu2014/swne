################################################################################
### Alignment workflow for the four human pancreatic islet datasets
################################################################################

library(Seurat)
library(Matrix)

# Read in all four input expression matrices
celseq.data <- read.table("pancreas_multi_celseq_expression_matrix.txt.gz")
celseq2.data <- read.table("pancreas_multi_celseq2_expression_matrix.txt.gz")
fluidigmc1.data <- read.table("pancreas_multi_fluidigmc1_expression_matrix.txt.gz")
smartseq2.data <- read.table("pancreas_multi_smartseq2_expression_matrix.txt.gz")

# Convert to sparse matrices for efficiency
celseq.data <- as(as.matrix(celseq.data), "dgCMatrix")
celseq2.data <- as(as.matrix(celseq2.data), "dgCMatrix")
fluidigmc1.data <- as(as.matrix(fluidigmc1.data), "dgCMatrix")
smartseq2.data <- as(as.matrix(smartseq2.data), "dgCMatrix")

# Create and setup Seurat objects for each dataset
celseq <- CreateSeuratObject(raw.data = celseq.data)
celseq <- FilterCells(celseq, subset.names = "nGene", low.thresholds = 1750)
celseq <- NormalizeData(celseq)
celseq <- FindVariableGenes(celseq, do.plot = F, display.progress = F)
celseq <- ScaleData(celseq)
celseq@meta.data$tech <- "celseq"

celseq2 <- CreateSeuratObject(raw.data = celseq2.data)
celseq2 <- FilterCells(celseq2, subset.names = "nGene", low.thresholds = 2500)
celseq2 <- NormalizeData(celseq2)
celseq2 <- FindVariableGenes(celseq2, do.plot = F, display.progress = F)
celseq2 <- ScaleData(celseq2)
celseq2@meta.data$tech <- "celseq2"

fluidigmc1 <- CreateSeuratObject(raw.data = fluidigmc1.data)
fluidigmc1 <- NormalizeData(fluidigmc1)
fluidigmc1 <- FindVariableGenes(fluidigmc1, do.plot = F, display.progress = F)
fluidigmc1 <- ScaleData(fluidigmc1)
fluidigmc1@meta.data$tech <- "fluidigmc1"

smartseq2 <- CreateSeuratObject(raw.data = smartseq2.data)
smartseq2 <- FilterCells(smartseq2, subset.names = "nGene", low.thresholds = 2500)
smartseq2 <- NormalizeData(smartseq2)
smartseq2 <- FindVariableGenes(smartseq2, do.plot = F, display.progress = F)
smartseq2 <- ScaleData(smartseq2)
smartseq2@meta.data$tech <- "smartseq2"

# Determine genes to use for CCA, must be highly variable in at least 2 datasets
ob.list <- list(celseq, celseq2, fluidigmc1, smartseq2)
genes.use <- c()
for (i in 1:length(ob.list)) {
  genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(ob.list)) {
  genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}

# Run multi-set CCA
pancreas.integrated <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 20)

# CC Selection
MetageneBicorPlot(pancreas.integrated, grouping.var = "tech", dims.eval = 1:20)

# Run rare non-overlapping filtering
pancreas.integrated <- CalcVarExpRatio(object = pancreas.integrated, reduction.type = "pca",
                                       grouping.var = "tech", dims.use = 1:20)
pancreas.integrated <- SubsetData(pancreas.integrated, subset.name = "var.ratio.pca",
                                           accept.low = 0.5)

# Alignment
pancreas.integrated <- AlignSubspace(pancreas.integrated,
                                     reduction.type = "cca",
                                     grouping.var = "tech",
                                     dims.align = 1:20)

# t-SNE and Clustering
pancreas.integrated <- FindClusters(pancreas.integrated, reduction.type = "cca.aligned",
                                    dims.use = 1:20, save.SNN = T, resolution = 0.4)
pancreas.integrated <- RunTSNE(pancreas.integrated,
                               reduction.use = "cca.aligned",
                               dims.use = 1:20)

saveRDS(pancreas.integrated, file = "pancreas_integrated_seurat.RObj")

# Visualization
TSNEPlot(pancreas.integrated, do.label = T)
TSNEPlot(pancreas.integrated, group.by = "tech")
