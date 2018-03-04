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

## Split data
set.seed(3523523)
cells.train <- sample(colnames(counts), size = ncol(counts)/2)
cells.test <- colnames(counts)[!colnames(counts) %in% cells.train]

## Embed training dataset
norm.counts.train <- ScaleCounts(counts[,cells.train], method = "ft")
nmf.res.train <- RunNMF(norm.counts.train[var.genes,], k = k, alpha = 0, init = "ica", n.cores = n.cores, loss = loss)
nmf.scores.train <- nmf.res.train$H

pc.train.res <- prcomp(t(norm.counts.train[var.genes,]), rank = 20)
pc.scores.train <- t(pc.train.res$x)

snn.train <- CalcSNN(pc.scores.train, k = 30, k.scale = 10)
swne.embedding.train <- EmbedSWNE(nmf.scores.train, snn.train, alpha.exp = 1, snn.exp = 1,
                                  n_pull = 4, dist.use = "IC")
swne.embedding.train$H.coords$name <- ""

pdf("pbmc3k_swne_plot_train_nolabel.pdf", height = 6, width = 6)
PlotSWNE(swne.embedding.train, alpha.plot = 0.4, sample.groups = cell.clusters[cells.train], do.label = T,
         label.size = 0, pt.size = 1.5, show.legend = F, seed = seed)
dev.off()

## Project test data onto training embedding
norm.counts.test <- ScaleCounts(counts[,cells.test], method = "ft")
nmf.scores.test <- ProjectSamples(norm.counts.test[var.genes,], nmf.res$W, loss = loss, n.cores = n.cores)

pc.scores.test <- t(t(norm.counts.test[var.genes,]) %*% pc.train.res$rotation)

snn.test <- ProjectSNN(pc.scores.test, pc.scores.train, k = 30, k.scale = 10)
sample.coords.test <- ProjectSWNE(swne.embedding.train, nmf.scores.test, SNN.test = snn.test, 
                                  alpha.exp = 1, snn.exp = 1, n_pull = 4)

swne.embedding.full <- swne.embedding.train
swne.embedding.full$sample.coords <- sample.coords.test

pdf("pbmc3k_swne_plot_test_nolabel.pdf", height = 6, width = 6)
PlotSWNE(swne.embedding.full, alpha.plot = 0.4, sample.groups = cell.clusters[cells.test], do.label = T,
         label.size = 0, pt.size = 1.5, show.legend = F, seed = seed)
dev.off()

