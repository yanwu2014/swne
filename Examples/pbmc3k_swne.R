library(Seurat)
library(swne)

## Read in Seurat object
se.obj <- readRDS("pbmc3k_seurat.Robj")

## Pull out counts matrix
counts <- as(se.obj@raw.data, "dgCMatrix")
counts <- counts[,colnames(counts) %in% se.obj@cell.names]

## Filter raw counts
counts <- FilterData(counts, min.samples.frac = 0.0025, min.nonzero.features = 200, trim = 0.001)
dim(counts)

## Normalize and scale raw counts
norm.counts <- ScaleCounts(counts, batch = NULL, method = "ft", adj.var = T)

## Pull out overdispersed genes as defined by Seurat
var.genes <- se.obj@var.genes; var.genes <- var.genes[var.genes %in% rownames(norm.counts)]
length(var.genes)

## Pull out cell clusters and defined by Seurat
cell.clusters <- se.obj@ident; names(cell.clusters) <- se.obj@cell.names;
levels(cell.clusters)

## Unguided NMF
loss <- "mse" ## Loss function
k.range <- seq(1,10,1) ## Range of factors to iterate over
n.cores <- 16 ## Number of cores to use
seed <- 32566 ## Set seed for cluster colors

## Identify optimal number of factors
n.comp.res <- FindNumFactors(norm.counts[var.genes,], k.range = k.range, n.cores = n.cores, do.plot = T, loss = loss)
n.comp.res$k

## Run NMF
k <- 8
nmf.res <- RunNMF(norm.counts[var.genes,], k = k, alpha = 0, init = "ica", n.cores = n.cores, loss = loss)
nmf.scores <- nmf.res$H

## Run SWNE embedding
snn.matrix <- se.obj@snn[colnames(nmf.scores), colnames(nmf.scores)]
swne.embedding <- EmbedSWNE(nmf.scores, snn.matrix, alpha.exp = 1.0, snn.exp = 1,
                            n_pull = 4, dist.use = "IC")

## Plot SWNE
PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = cell.clusters, do.label = T,
         label.size = 4, pt.size = 1.5, show.legend = T, seed = seed)

## Plot tSNE for comparison
tsne.scores <- GetCellEmbeddings(se.obj, reduction.type = "tsne")
PlotDims(tsne.scores, sample.groups = cell.clusters, pt.size = 1, label.size = 3, alpha = 0.5,
         show.legend = F, seed = seed)

## Annotate factors with marker genes
gene.factor.cors <- FactorAssociation(norm.counts, nmf.scores, n.cores = n.cores, metric = "pearson")
top.factor.genes.df <- SummarizeAssocFeatures(gene.factor.cors, features.return = 5)

## Plot heatmap
gene.factor.cors.plot <- gene.factor.cors[unique(top.factor.genes.df$feature),]
ggHeat(gene.factor.cors.plot, clustering = "col")

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
         label.size = 4, pt.size = 1.5, show.legend = T, seed = seed)
