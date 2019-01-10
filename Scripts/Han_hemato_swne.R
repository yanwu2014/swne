library(Matrix)
library(swne)

## Load data
load("Han/BM_BMcKit_PB_RData/xp.RData")
load("Han/BM_BMcKit_PB_RData/g2.RData")
load("Han/BM_BMcKit_PB_RData/cells_AUC.RData")
hemato.genes <- read.table("Han/annots/haematopedia_mmc6.csv", header = T, sep = ",")

## Filter dataset
w2 <- !is.na(g2)
xp <- Matrix::t(xp[w2,])

## Assign labels to cells
lineages <- c("Multi Potential Progenitor", "Macrophage Lineage", "Neutrophil Lineage",
              "Erythrocyte Lineage", "B Cell Lineage", "T Cell Lineage", "NK Cell Lineage")
cutoffs <- setNames(c(0.04,0.09,0.05,0.045,0.09,0.075,0.04), lineages)

labels <- sapply(lineages, function(i) cells_AUC@assays[[1]][i,][w2] >= cutoffs[i])
labels <- apply(labels, 1, which)
labels <- sapply(labels, function(x) { if(length(x) == 1) {x} else {0} })
labels[labels != 0] <- lineages[labels[labels != 0]]
labels[labels == 0] <- NA
names(labels) <- colnames(xp)
labels <- factor(labels)
labels <- plyr::revalue(labels, replace = c("Multi Potential Progenitor" = "MPP",
                                            "Macrophage Lineage" = "Macrophage",
                                            "Neutrophil Lineage" = "Neutrophil",
                                            "Erythrocyte Lineage" = "Erythrocyte",
                                            "B Cell Lineage" = "B Cell",
                                            "T Cell Lineage" = "T Cell",
                                            "NK Cell Lineage" = "NK Cell"))
table(labels); paste("Cells with missing labels:", sum(is.na(labels)));

## Run PCA
load("Han/BM_BMcKit_PB_RData/pca_g2.RData")
rownames(pca) <- names(labels)

## Make t-SNE and UMAP plots
plot.seed <- 312525

## Run t-SNE
load("Han/BM_BMcKit_PB_RData/tsne_g2.RData")
rownames(tsne) <- names(labels)
pdf("Han_hemato_tsne_plot.pdf", width = 6, height = 6)
PlotDims(tsne, sample.groups = labels, show.legend = F, show.axes = F,
         alpha.plot = 0.75, label.size = 4, pt.size = 0.5,
         seed = plot.seed, use.brewer.pal = T)
dev.off()

## Run UMAP
load("Han/BM_BMcKit_PB_RData/umap_g2.RData")
rownames(umap) <- names(labels)
pdf("Han_hemato_umap_plot.pdf", width = 6, height = 6)
PlotDims(umap, sample.groups = labels, show.legend = F, show.axes = F,
         alpha.plot = 0.75, label.size = 4, pt.size = 0.5,
         seed = plot.seed, use.brewer.pal = T)
dev.off()

## Run SWNE

## Filter lowly expressed genes and get gene variance info
var.df <- AdjustVariance(xp, verbose = F, plot = T)

n.genes <- 2e3
var.df <- var.df[order(var.df$lp),]
var.genes <- rownames(var.df[1:n.genes,])

n.cores <- 24
nmf.res <- RunNMF(xp[var.genes,], k = 20, n.cores = n.cores, ica.fast = T)
nmf.res$W <- ProjectFeatures(xp, nmf.res$H, n.cores = n.cores)

snn <- CalcSNN(t(pca), k = 40, prune.SNN = 0.0)
swne.embedding <- EmbedSWNE(nmf.res$H, SNN = snn, alpha.exp = 1.5, snn.exp = 0.5, n_pull = 3)
swne.embedding$H.coords$name <- ""

## Embed selected genes onto swne plot
genes.embed <- c("Ms4a1", "Cd4", "Ly6g", "Nkg7", "Fcgr1")
swne.embedding <- EmbedFeatures(swne.embedding, nmf.res$W, genes.embed,
                                n_pull = 3)

## SWNE plots
pdf("Han_hemato_swne_plot.pdf", width = 6, height = 6)
PlotSWNE(swne.embedding, alpha.plot = 0.6, sample.groups = labels, do.label = T,
         label.size = 4, pt.size = 0.75, show.legend = F, seed = plot.seed,
         use.brewer.pal = T)
dev.off()


## Quantitative evaluation of t-SNE, UMAP, SWNE
library(FNN)
library(proxy)

## Calculate approximate kNN for an embedding
ComputeKNN <- function(emb, k) {
  knn.idx <- knn.index(t(emb), k = k)
  knn.matrix <- matrix(0, ncol(emb), ncol(emb))
  for (i in 1:nrow(knn.idx)) {
    knn.matrix[knn.idx[i,],i] <- 1
    knn.matrix[i, knn.idx[i,]] <- 1
  }
  rownames(knn.matrix) <- colnames(knn.matrix) <- colnames(emb)
  as(knn.matrix, "dgCMatrix")
}


## Calculate Jaccard similarities
CalcJaccard <- function(x,y) {
  a <- sum(x)
  b <- sum(y)
  c <- sum(x == 1 & y == 1)
  c/(a + b - c)
}


## Calculate pairwise distances between centroids
CalcPairwiseDist <- function(data.use, clusters, dist.method = "euclidean") {
  data.centroids <- t(apply(data.use, 1, function(x) tapply(x, clusters, mean)))
  return(proxy::dist(data.centroids, method = dist.method, by_rows = F))
}


## Compile embeddings
embeddings <- list(swne = t(as.matrix(swne.embedding$sample.coords)),
                   tsne = t(tsne), umap = t(umap))


## Compute cluster distance correlations
label.cells <- names(labels[!is.na(labels)])
ref.dist <- CalcPairwiseDist(xp[,label.cells], labels[label.cells])

embeddings.cor <- sapply(embeddings, function(emb) {
  emb.dist <- CalcPairwiseDist(emb[,label.cells], labels[label.cells])
  cor(ref.dist, emb.dist)
})
print(embeddings.cor)

## Calculate neighborhood fidelity
n.neighbors <- 30
ref.knn <- ComputeKNN(xp[,label.cells], k = n.neighbors)

## Compute kNN for embeddings
embeddings.knn <- lapply(embeddings, function(x) ComputeKNN(x[,label.cells], k = n.neighbors))

knn.simil <- sapply(embeddings.knn, function(knn.emb) {
  mean(sapply(1:ncol(knn.emb), function(i) CalcJaccard(knn.emb[,i], ref.knn[,i])))
})
print(knn.simil)


library(ggplot2)
library(ggrepel)

pdf("Han_hemato_quant_eval.pdf", width = 5, height = 4.5)
scatter.df <- data.frame(x = knn.simil, y = embeddings.cor, name = names(embeddings))
ggplot(scatter.df, aes(x, y)) + geom_point(size = 2, alpha = 1) +
  theme_classic() + theme(legend.position = "none", text = element_text(size = 16)) +
  xlab("Neighborhood Similarity") + ylab("Cluster Distance Correlation") +
  geom_text_repel(aes(x, y, label = name), size = 6.5) +
  xlim(0, max(knn.simil)) + ylim(0, max(embeddings.cor))
dev.off()

save.image("Han_hemato_swne.RData")
