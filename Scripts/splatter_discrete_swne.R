#### Run SWNE ####
library(swne)

## Set directory
setwd("/media/Scratch_SSD/Yan/R_Analysis/swne-dev/simulations")

## Load data
load("splatter.discrete.RData")

## Pull out relevant variables
pc.scores <- r$reductions$PCA
tsne.scores <- r$embeddings$PCA$tSNE
var.genes <- rownames(r$misc$PCA$v)

## Pull out group info
clusters <- factor(metadata$Group); names(clusters) <- rownames(r$misc$rawCounts);

## Run SWNE
loss <- "mse"
nmf.init <- "ica"
n.cores <- 32

norm.counts <- ScaleCounts(Matrix::t(r$misc$rawCounts), batch = NULL, method = "log", adj.var = T)
snn <- CalcSNN(t(pc.scores), k = 30, prune.SNN = 0.0)

k <- 6
nmf.res <- RunNMF(norm.counts[var.genes,], k = k, alpha = 0, init = nmf.init, 
                  n.cores = n.cores, loss = loss)
H <- nmf.res$H

nmf.res$W <- ProjectFeatures(norm.counts, nmf.res$H, loss = "mse", n.cores = n.cores)
top.genes.df <- SummarizeAssocFeatures(nmf.res$W, features.return = 2)

swne.embedding <- EmbedSWNE(H, SNN = snn, alpha.exp = 1.0, snn.exp = 0.25, n_pull = 4,
                            snn.factor.proj = T, dist.use = "cosine")
swne.embedding <- EmbedFeatures(swne.embedding, nmf.res$W, c("Gene2675", "Gene25475"), n_pull = k)
swne.embedding$H.coords$name <- ""

swne.embedding.no.snn <- EmbedSWNE(H, SNN = NULL, alpha.exp = 1.0, snn.exp = 0.25, n_pull = 4,
                                   dist.use = "cosine")
swne.embedding.no.snn <- EmbedFeatures(swne.embedding.no.snn, nmf.res$W, c("Gene2675", "Gene25475"), n_pull = k)
swne.embedding.no.snn$H.coords$name <- ""

#### Generate plots ####

## Make SWNE plot
color.seed <- 43859279
pdf("splatter_discrete_swne.pdf", width = 4, height = 4)
PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = clusters, do.label = T, label.size = 4, 
         pt.size = 1.5, seed = color.seed, show.legend = F)
dev.off()

pdf("splatter_discrete_swne_no_snn.pdf", width = 4, height = 4)
PlotSWNE(swne.embedding.no.snn, alpha.plot = 0.4, sample.groups = clusters, do.label = T, label.size = 4, 
         pt.size = 1.5, seed = color.seed, show.legend = F)
dev.off()

## Make PCA plots
pdf("splatter_discrete_pca1-2.pdf", width = 4, height = 4)
PlotDims(pc.scores[,1:2], sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()

pdf("splatter_discrete_pca2-3.pdf", width = 4, height = 4)
PlotDims(pc.scores[,c(2,3)], sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()


## Make t-SNE plot
pdf("splatter_discrete_tsne.pdf", width = 4, height = 4)
PlotDims(tsne.scores, sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()


## Make MDS plot (calculate distance matrix first)
mds.scores <- cmdscale(dist(t(as.matrix(norm.counts[var.genes,]))), k = 2)
pdf("splatter_discrete_mds.pdf", width = 4, height = 4)
PlotDims(mds.scores, sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 3.5, show.legend = F)
dev.off()


## Make Isomap plot
library(RDRToolbox)
isomap.scores <- Isomap(t(as.matrix(norm.counts[var.genes,])), dims = 2, k = 20)$dim2
pdf("splatter_discrete_isomap.pdf", width = 4, height = 4)
PlotDims(isomap.scores, sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()


## Make LLE plot
lle.scores <- LLE(t(as.matrix(norm.counts[var.genes,])), dim = 2, k = 20)
pdf("splatter_discrete_lle.pdf", width = 4, height = 4)
PlotDims(lle.scores, sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()


## Make diffusion maps
library(destiny)
dm <- DiffusionMap(t(as.matrix(norm.counts[var.genes,])), k = 20, n_eigs = 8)
diff.scores <- dm@eigenvectors

pdf("splatter_discrete_dmaps.pdf", width = 4, height = 4)
PlotDims(diff.scores[,c(1,2)], sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()


## Make UMAP plot
# write.table(as.matrix(norm.counts[var.genes,]), file = "splatter.discrete.norm.counts.tsv", sep = "\t")
umap.scores <- as.matrix(read.table("splatter.discrete.umap.tsv", sep = "\t", header = T, row.names = 1))
rownames(umap.scores) <- colnames(norm.counts)

pdf("splatter_discrete_umap_genes.pdf", width = 5, height = 5)
PlotDims(umap.scores, sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.3, label.size = 3.5, show.legend = F)
dev.off()

save.image("splatter_discrete_analysis.RData")


#### Evaluate embeddings ####

## Set directory
setwd("/media/Scratch_SSD/Yan/R_Analysis/swne-dev/simulations")
load("splatter_discrete_analysis.RData")
library(swne)
library(ggplot2)
source("splatter_helper_functions.R")

## Calculate pairwise distances in original space
clust.dist <- CalcPairwiseDist(norm.counts[var.genes,], clusters, "euclidean")

## Correlate pairwise distances for embeddings
embeddings <- list(swne = t(as.matrix(swne.embedding$sample.coords)),
                   swne.no.snn = t(as.matrix(swne.embedding.no.snn$sample.coords)),
                   tsne = t(tsne.scores),
                   pca = t(pc.scores[,1:2]),
                   mds = t(mds.scores),
                   lle = t(lle.scores),
                   isomap = t(isomap.scores),
                   dmaps = t(diff.scores))

embeddings.cor <- sapply(embeddings, function(emb) {
  emb.dist <- CalcPairwiseDist(emb, clusters, "euclidean")
  cor(emb.dist, clust.dist)
})
print(embeddings.cor)


## Calculate silhouette scores for embeddings
si.avg.scores <- sapply(embeddings, function(emb) {
  emb.dist <- dist(t(emb))
  mean(cluster::silhouette(as.integer(clusters), emb.dist)[,3])
})
print(si.avg.scores)


## Plot local and global evaluations together
library(ggplot2)
library(ggrepel)
scatter.df <- data.frame(x = si.avg.scores, y = embeddings.cor, name = names(embeddings))

pdf("splatter_discrete_quant_eval.pdf", width = 4.5, height = 4)
ggplot(scatter.df, aes(x, y)) + geom_point(size = 2, alpha = 1) +
  theme_classic() + theme(legend.position = "none", text = element_text(size = 14)) +
  xlab("Silhouette Score") + ylab("Cluster Distance Correlation") + 
  geom_text_repel(aes(x, y, label = name), size = 4.5) + 
  xlim(0, max(si.avg.scores)) + ylim(0, max(embeddings.cor))
dev.off()

write.table(scatter.df, file = "splatter_discrete_quant_eval.tsv", sep = "\t")

save.image("splatter_discrete_analysis.RData")


#### Evaluate how SWNE visualization changes with k ####

setwd("/media/Scratch_SSD/Yan/R_Analysis/swne-dev/simulations")
load("splatter_discrete_analysis.RData")
source("splatter_helper_functions.R")
library(swne)

k.range <- c(3, 4, 6, 8, 10, 12)
k.range.embeddings <- lapply(k.range, function(k) {
  nmf.res <- RunNMF(norm.counts[var.genes,], k = k, alpha = 0, init = "ica", 
                    n.cores = n.cores, loss = loss)
  swne.embedding <- EmbedSWNE(nmf.res$H, SNN = snn, alpha.exp = 1.0, snn.exp = 0.25, n_pull = 4,
                              snn.factor.proj = T, dist.use = "cosine")
  swne.embedding$H.coords$name <- ""
  
  pdf(paste("splatter_discrete_swne_k", k, ".pdf", sep = ""), width = 3.5, height = 3.5)
  print(PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = clusters, do.label = T, label.size = 4.5,
                 pt.size = 1.5, seed = color.seed, show.legend = F))
  dev.off()
  
  return(t(as.matrix(swne.embedding$sample.coords)))
})
names(k.range.embeddings) <- paste0("k = ", k.range)


## Correlate pairwise distances for embeddings
k.range.embeddings.cor <- sapply(k.range.embeddings, function(emb) {
  emb.dist <- CalcPairwiseDist(emb, clusters, "euclidean")
  cor(emb.dist, clust.dist)
})
print(k.range.embeddings.cor)

## Calculate silhouette scores for embeddings
k.range.si.avg.scores <- sapply(k.range.embeddings, function(emb) {
  emb.dist <- dist(t(emb))
  mean(cluster::silhouette(as.integer(clusters), emb.dist)[,3])
})
print(k.range.si.avg.scores)


pdf("k_range_discrete_cl_dist_cor.pdf", width = 3, height = 2.5)
ggLinePlot(k.range.embeddings.cor)
dev.off()

pdf("k_range_discrete_cl_silhouette.pdf", width = 3, height = 2.5)
ggLinePlot(k.range.si.avg.scores)
dev.off()

save.image("splatter_discrete_analysis.RData")
