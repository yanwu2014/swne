library(pagoda2)
library(swne)
library(Matrix)

## Load pagoda2 object
p2 <- readRDS("pbmc3k_pagoda2.Robj")

## Pull out raw counts and filter
norm.counts <- ScaleCounts(t(p2$misc$rawCounts[ ,colnames(p2$counts)]), method = "log", adj.var = T)
dim(norm.counts)

## Pull out clusters
clusters <- p2$clusters$PCA$multilevel; clusters <- clusters[colnames(norm.counts)]

## Pull out variable genes
n.od.genes <- 3e3
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
od.genes <- rownames(var.info[1:n.od.genes,])
length(od.genes)

## NMF parameters
loss <- "mse"
n.cores <- p2$n.cores
init <- "ica" # Initialization method for NMF

# ## Find number of factors to use
# k.res <- FindNumFactors(norm.counts[od.genes,], k.range = seq(2,12,2), n.cores = n.cores, do.plot = T,
#                         na.frac = 0.25)

## Run NMF
k <- 10
nmf.res <- RunNMF(norm.counts[od.genes, ], k = k, init = init, n.cores = n.cores, loss = loss,
                  init.zeros = "random")
nmf.res$W <- ProjectFeatures(norm.counts, nmf.res$H, loss = "mse", n.cores = n.cores)
nmf.scores <- nmf.res$H ## Pull out NMF components

## SWNE embedding parameters
alpha.exp <- 1.5 # Increase this to move the cells closer to the factors. Values > 2 start to distort the data.
snn.exp <- 0.5 # Lower this to move similar cells closer to each other
n_pull <- 4 # The number of factors pulling on each cell. Must be at least 3.
dist.use <- "IC" # You can also use pearson correlation, but IC seems to work better.

## Build SNN matrix from PC scores
pc.scores <- t(p2$reductions$PCA)
snn <- CalcSNN(pc.scores, k = 30, k.scale = 10, prune.SNN = 1/15)

## Run SWNE embedding
swne.embedding <- EmbedSWNE(nmf.scores, snn, alpha.exp = alpha.exp, snn.exp = snn.exp, n_pull = n_pull,
                            dist.use = dist.use, min.snn = 0.0)

## Set seed for color shuffling
color.seed <- 32566

## SWNE plot
PlotSWNE(swne.embedding, alpha.plot = 0.3, sample.groups = clusters, do.label = T, label.size = 3.5,
         pt.size = 1.25, show.legend = F, seed = color.seed)

## Plot tSNE for comparison
tsne.scores <- p2$embeddings$PCA$tSNE
PlotDims(tsne.scores, sample.groups = clusters, pt.size = 0.75, label.size = 3, alpha = 0.3,
         show.legend = F, seed = color.seed, show.axes = F)

## Hide factors for now
swne.embedding$H.coords$name <- ""

## Embed key features
genes.embed <- c("MS4A1", "GNLY", "CD3E", "CD14",
                 "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
swne.embedding <- EmbedFeatures(swne.embedding, nmf.res$W, genes.embed,
                                alpha.exp = 1, n_pull = 4, scale.cols = T)

## SWNE plot
PlotSWNE(swne.embedding, alpha.plot = 0.3, sample.groups = clusters, do.label = T, label.size = 3.5,
         pt.size = 1.25, show.legend = F, seed = color.seed)

## Annotate factors with marker genes
gene.loadings <- nmf.res$W
top.factor.genes.df <- SummarizeAssocFeatures(gene.loadings, features.return = 2)

## Plot heatmap
gene.loadings.heat <- gene.loadings[unique(top.factor.genes.df$feature),]
ggHeat(gene.loadings.heat, clustering = "both")

