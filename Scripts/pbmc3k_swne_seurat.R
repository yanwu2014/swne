library(Seurat)
library(swne)

## Read in Seurat object
se.obj <- readRDS("pbmc3k_seurat.Robj")

## Pull out counts matrix, filter and trim data
counts <- as(se.obj@raw.data, "dgCMatrix")[,se.obj@cell.names]
counts <- FilterData(counts, min.samples.frac = 0.0025, min.nonzero.features = 200, trim = 0.001)
dim(counts)

## Normalize and scale raw counts
norm.counts <- ScaleCounts(counts, batch = NULL, method = "log", adj.var = T)

## Pull out overdispersed genes as defined by Seurat
var.genes <- intersect(se.obj@var.genes, rownames(counts));
length(var.genes)

## Pull out cell clusters as defined by Seurat
cell.clusters <- se.obj@ident; names(cell.clusters) <- se.obj@cell.names;
levels(cell.clusters)

## Unguided NMF
loss <- "mse" ## Loss function
n.cores <- 16 ## Number of cores to use

## Identify optimal number of factors
k.range <- seq(1,10,1) ## Range of factors to iterate over
k.res <- FindNumFactors(norm.counts[var.genes,], k.range = k.range, n.cores = n.cores, do.plot = T, loss = loss)
k.res$k

## Run NMF
k <- 8
nmf.res <- RunNMF(norm.counts[var.genes,], k = k, alpha = 0, init = "ica", n.cores = n.cores, loss = loss)
nmf.scores <- nmf.res$H

## Compute snn or use pre-computed snn from Seurat
# pc.scores <- t(GetCellEmbeddings(se.obj, reduction.type = "pca", dims.use = 1:k))
# snn <- CalcSNN(pc.scores)
snn <- se.obj@snn

## Run SWNE embedding
alpha.exp <- 1.5 # Increase this to move the cells closer to the factors. Values > 2 start to distort the data.
snn.exp <- 0.5 # Lower this to move similar cells closer to each other
n_pull <- 4 # The number of factors pulling on each cell. Must be at least 3.
swne.embedding <- EmbedSWNE(nmf.scores, snn, alpha.exp = alpha.exp, snn.exp = snn.exp,
                            n_pull = n_pull, dist.use = "IC")

## Plot SWNE
color.seed <- 32566 ## Set seed for cluster colors
PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = cell.clusters, do.label = T,
         label.size = 3.5, pt.size = 1.25, show.legend = T, seed = color.seed)

## Plot tSNE for comparison
tsne.scores <- GetCellEmbeddings(se.obj, reduction.type = "tsne")
PlotDims(tsne.scores, sample.groups = cell.clusters, pt.size = 1, label.size = 3.5, alpha = 0.4,
         show.legend = F, seed = color.seed, show.axes = F)

## Annotate factors with marker genes
gene.loadings <- nmf.res$W
top.factor.genes.df <- SummarizeAssocFeatures(gene.loadings, features.return = 2)

## Plot heatmap
gene.loadings.heat <- gene.loadings[unique(top.factor.genes.df$feature),]
ggHeat(gene.loadings.heat, clustering = "col")

## Rename some factors
# swne.embedding <- RenameFactors(swne.embedding, name.mapping = c("factor_1" = "B-Cell Signaling"))
swne.embedding$H.coords$name <- ""

## Pick some key genes to embed
genes.embed <- c("MS4A1", "GNLY", "CD3E", "CD14",
                 "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")

## Calculate gene loadings for all genes
nmf.res$W <- ProjectFeatures(norm.counts, nmf.scores, loss = loss, n.cores = n.cores)

## Project key genes onto embedding
swne.embedding <- EmbedFeatures(swne.embedding, nmf.res$W, genes.embed, n_pull = 3)

## Remake SWNE plot
PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = cell.clusters, do.label = T,
         label.size = 4, pt.size = 1.5, show.legend = T, seed = color.seed)
