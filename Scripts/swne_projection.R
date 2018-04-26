## Demonstrates how to project new data onto an existing SWNE embedding
library(Seurat)
library(swne)

#### Split pbmc3k dataset in half and project one half onto the other ####

## Read in Seurat object
se.obj <- readRDS("pbmc3k_seurat.Robj")

## Pull out counts matrix
counts <- FilterData(as(se.obj@raw.data, "dgCMatrix")[,se.obj@cell.names], min.samples.frac = 0.005,
                     min.nonzero.features = 200, trim = 5)
dim(counts)

## Split data
set.seed(3523523)
cells.train <- sample(colnames(counts), size = ncol(counts)/2)
cells.test <- colnames(counts)[!colnames(counts) %in% cells.train]

norm.counts.train <- ScaleCounts(counts[,cells.train], method = "ft")

## Embed training dataset
k <- 10
loss <- "mse"
n.cores <- 16
var.genes <- intersect(se.obj@var.genes, rownames(counts))

nmf.res <- RunNMF(norm.counts.train[var.genes,], k = k, init = "ica", n.cores = n.cores, loss = loss)

pc.train <- prcomp(t(norm.counts.train[var.genes,]), center = T, scale = F, rank = 20)
pc.scores <- t(pc.train$x)

snn <- CalcSNN(pc.scores, k = 30)
swne.embedding.train <- EmbedSWNE(nmf.res$H, snn, alpha.exp = 1.5, snn.exp = 0.25,
                                  n_pull = 4, dist.use = "cosine")
swne.embedding.train$H.coords$name <- ""

clusters <- se.obj@ident; names(clusters) <- se.obj@cell.names;
color.seed <- 3362083

pdf("pbmc3k_swne_plot_train.pdf", height = 6, width = 6)
PlotSWNE(swne.embedding.train, alpha.plot = 0.4, sample.groups = clusters, do.label = T,
         label.size = 4.5, pt.size = 1.25, show.legend = F, seed = color.seed)
dev.off()

## Project test data onto training embedding
norm.counts.test <- ScaleCounts(counts[,cells.test], method = "ft")
H.test <- ProjectSamples(norm.counts.test[var.genes,], nmf.res$W, loss = loss, n.cores = n.cores)

test.cv <- Matrix::rowMeans(norm.counts.test[var.genes,])
pc.scores.test <- t((t(norm.counts.test[var.genes,]) - test.cv) %*% pc.train$rotation)

snn.test <- ProjectSNN(pc.scores.test, pc.scores, k = 20, k.scale = 10)
swne.embedding.test <- ProjectSWNE(swne.embedding.train, H.test, SNN = snn.test, alpha.exp = 1.5,
                                   snn.exp = 0.25, n_pull = 4)

pdf("pbmc3k_swne_plot_test.pdf", height = 6, width = 6)
PlotSWNE(swne.embedding.test, alpha.plot = 0.4, sample.groups = clusters, do.label = T,
         label.size = 4.5, pt.size = 1.25, show.legend = F, seed = color.seed)
dev.off()


#### Project pbmc33k dataset onto pbmc3k dataset ####

## Project pbmc33k dataset
pbmc33k <- readRDS("pbmc33k_seurat.Robj")
norm.counts.33k <- ScaleCounts(as(pbmc33k@raw.data[,pbmc33k@cell.names], "dgCMatrix"), method = "ft")

## Project NMF
proj.genes.33k <- intersect(var.genes, rownames(norm.counts.train))
H.33k <- ProjectSamples(norm.counts.33k[proj.genes.33k,], nmf.res$W[proj.genes.33k,], loss = loss, n.cores = n.cores)

## Project SNN
cm <- rowMeans(norm.counts.33k[proj.genes.33k,])
pc.scores.33k <- t(t(norm.counts.33k[proj.genes.33k,] - cm) %*% pc.train$rotation[proj.genes.33k,])

snn.33k <- ProjectSNN(pc.scores.33k, pc.scores, k = 20, k.scale = 10)

## Create test SWNE embedding
swne.embedding.33k <- ProjectSWNE(swne.embedding.train, H.33k, snn.33k, alpha.exp = 1.5,
                                   snn.exp = 0.25, n_pull = 4)

clusters.33k <- pbmc33k@ident; names(clusters.33k) <- pbmc33k@cell.names;

pdf("pbmc33k_swne_plot_projected.pdf", height = 7, width = 8)
PlotSWNE(swne.embedding.33k, alpha.plot = 0.3, sample.groups = clusters.33k, do.label = T,
         label.size = 4, pt.size = 1.0, show.legend = F, seed = color.seed)
dev.off()

save.image("swne_projection_demo.RData")
