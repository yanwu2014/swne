#### Run SWNE ####

library(swne)

## Set directory
setwd("/media/Scratch_SSD/Yan/R_Analysis/swne-dev/simulations")

## Load data
load("splatter.trajectory.RData")

## Pull out relevant variables
pc.scores <- r$reductions$PCA
tsne.scores <- r$embeddings$PCA$tSNE
clusters <- metadata$Group; names(clusters) <- colnames(counts);

## Pull out variable genes
var.genes <- rownames(r$misc$PCA$v)

## Run SWNE
loss <- "mse"
nmf.init <- "random"
n.cores <- 32

## Scale counts
norm.counts <- ScaleCounts(Matrix::t(r$misc$rawCounts), batch = NULL, method = "log", adj.var = T)

## Make SNN
snn <- CalcSNN(t(pc.scores), k = 30, k.scale = 10, prune.SNN = 0.05)

## Run SWNE
nmf.res <- RunNMF(norm.counts[var.genes,], k = 10, alpha = 0, init = "random", 
                  n.cores = n.cores, loss = loss, n.rand.init = 5)
H <- nmf.res$H

nmf.res$W <- ProjectFeatures(norm.counts, nmf.res$H, loss = loss, n.cores = n.cores)
top.genes.df <- SummarizeAssocFeatures(nmf.res$W, features.return = 2)

swne.embedding <- EmbedSWNE(H, SNN = snn, alpha.exp = 2.5, snn.exp = 0.1, n_pull = 4, snn.factor.proj = T, 
                            pca.red = T, dist.use = "cosine")
swne.embedding <- EmbedFeatures(swne.embedding, nmf.res$W, c("Gene418", "Gene6053"), n_pull = 4)
swne.embedding$H.coords$name <- ""

swne.embedding.no.snn <- EmbedSWNE(H, SNN = NULL, alpha.exp = 2.5, snn.exp = 0.1, n_pull = 4,
                                   pca.red = T, dist.use = "cosine")
swne.embedding.no.snn <- EmbedFeatures(swne.embedding.no.snn, nmf.res$W, c("Gene418", "Gene6053"), n_pull = 4)
swne.embedding.no.snn$H.coords$name <- ""

#### Generate plots ####

## Make SWNE plot
color.seed <- 43859279
pdf("splatter_trajectory_swne.pdf", width = 4, height = 4)
PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = clusters, do.label = T, label.size = 4, 
         pt.size = 1.5, seed = color.seed, show.legend = F)
dev.off()

pdf("splatter_trajectory_swne_no_snn.pdf", width = 4, height = 4)
PlotSWNE(swne.embedding.no.snn, alpha.plot = 0.4, sample.groups = clusters, do.label = T, label.size = 4, 
         pt.size = 1.5, seed = color.seed, show.legend = F)
dev.off()


## Make PCA plots
pdf("splatter_trajectory_pca1-2.pdf", width = 4, height = 4)
PlotDims(pc.scores[,c(1,2)], sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()

pdf("splatter_trajectory_pca1-3.pdf", width = 4, height = 4)
PlotDims(pc.scores[,c(1,3)], sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()


## Make t-SNE plot
pdf("splatter_trajectory_tsne.pdf", width = 4, height = 4)
PlotDims(tsne.scores, sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()


## Make MDS plot (calculate distance matrix first)
mds.scores <- cmdscale(dist(t(as.matrix(norm.counts[var.genes,]))), k = 2)
pdf("splatter_trajectory_mds.pdf", width = 4, height = 4)
PlotDims(mds.scores, sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()


## Make Isomap plot
library(RDRToolbox)
isomap.scores <- Isomap(t(as.matrix(norm.counts[var.genes,])), dims = 2, k = 20)$dim2
rownames(isomap.scores) <- colnames(norm.counts)

pdf("splatter_trajectory_isomap.pdf", width = 4, height = 4)
PlotDims(isomap.scores, sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()


## Make LLE plot
lle.scores <- LLE(t(as.matrix(norm.counts[var.genes,])), dim = 2, k = 20)
rownames(lle.scores) <- colnames(norm.counts)

pdf("splatter_trajectory_lle.pdf", width = 4, height = 4)
PlotDims(lle.scores, sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()


## Make diffusion maps
library(destiny)
dm <- DiffusionMap(t(as.matrix(norm.counts[var.genes,])), k = 20, n_eigs = 6)
diff.scores <- dm@eigenvectors
rownames(diff.scores) <- colnames(norm.counts)

pdf("splatter_trajectory_dmaps.pdf", width = 4, height = 4)
PlotDims(diff.scores[,c(1,2)], sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()


## Make UMAP plot using genes
# write.table(as.matrix(norm.counts[var.genes,]), file = "splatter.trajectory.norm.counts.tsv", sep = "\t")
umap.scores <- as.matrix(read.table("splatter.trajectory.umap.tsv", sep = "\t", header = T, row.names = 1))
rownames(umap.scores) <- colnames(norm.counts)

pdf("splatter_trajectory_umap_genes.pdf", width = 4, height = 4)
PlotDims(umap.scores, sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()


## Make UMAP plot using PCs
# write.table(t(pc.scores), file = "splatter.trajectory.pcs.tsv", sep = "\t")
umap.scores <- as.matrix(read.table("splatter.trajectory.umap.pcs.tsv", sep = "\t", header = T, row.names = 1))
rownames(umap.scores) <- colnames(norm.counts)

pdf("splatter_trajectory_umap_pcs.pdf", width = 4, height = 4)
PlotDims(umap.scores, sample.groups = clusters, x.lab = NULL, y.lab = NULL, seed = color.seed, 
         show.axes = T, alpha.plot = 0.6, label.size = 4, show.legend = F)
dev.off()

save.image("splatter_trajectory_analysis.RData")


#### Evaluate embeddings ####

setwd("/media/Scratch_SSD/Yan/R_Analysis/swne-dev/simulations")
load("splatter_trajectory_analysis.RData")
library(swne)
library(ggplot2)
source("splatter_helper_functions.R")

## Make a list of all embeddings
embeddings <- list(swne = t(as.matrix(swne.embedding$sample.coords)),
                   swne.no.snn = t(as.matrix(swne.embedding.no.snn$sample.coords)),
                   tsne = t(tsne.scores),
                   pca = t(pc.scores[,1:2]),
                   mds = t(mds.scores),
                   lle = t(lle.scores),
                   isomap = t(isomap.scores),
                   dmaps = t(diff.scores[,1:2]))

## Create reference graph for trajectory simulation
traj.nn <- matrix(0, nrow(metadata), nrow(metadata))
rownames(traj.nn) <- colnames(traj.nn) <- rownames(metadata)

path.step <- GetPathStep(metadata, step.size = 1, make.factor = F)

## Connect cells within paths
for(i in 1:ncol(traj.nn)) {
  ps <- path.step[[i]]
  ps.cells <- names(path.step[path.step == ps])
  traj.nn[ps.cells, i] <- 1; traj.nn[i, ps.cells] <- 1;
  
  p <- strsplit(ps, split = "\\.")[[1]][[1]]
  s <- as.integer(strsplit(ps, split = "\\.")[[1]][[2]])
  if (s > 0) {
    ps.nn <- paste(p, s - 1, sep = ".")
    ps.nn.cells <- names(path.step[path.step == ps.nn])
    traj.nn[ps.nn.cells, i] <- 1; traj.nn[i, ps.nn.cells] <- 1;
    if (s > 1 && length(ps.nn.cells) < 5) {
      ps.nn2 <- paste(p, s - 2, sep = ".")
      ps.nn2.cells <- names(path.step[path.step == ps.nn2])
      traj.nn[ps.nn2.cells, i] <- 1; traj.nn[i, ps.nn2.cells] <- 1;
    }
  }
}

## Connect different paths
## Connect the end of path 1 with the start of path 2
path1.end.cells <- names(path.step[path.step == "Path1.100"])
path2.start.cells <- names(path.step[path.step == "Path2.1"])
traj.nn[path1.end.cells, path2.start.cells] <- 1; traj.nn[path2.start.cells, path1.end.cells] <- 1;

## Connect the end of path 1 with the start of path 3
path3.start.cells <- names(path.step[path.step == "Path3.1"])
traj.nn[path1.end.cells, path3.start.cells] <- 1; traj.nn[path3.start.cells, path1.end.cells] <- 1;

## Connect the end of path 3 with the start of path 4
path3.end.cells <- names(path.step[path.step == "Path3.50"])
path4.start.cells <- names(path.step[path.step == "Path4.1"])
traj.nn[path3.end.cells, path4.start.cells] <- 1; traj.nn[path4.start.cells, path3.end.cells] <- 1;

## Convert to sparse matrix to save space
traj.nn <- as(traj.nn, "dgCMatrix")

# ## Visualize reference NN graph to ensure it's correct
# library(igraph)
# library(RColorBrewer)
# 
# traj.graph <- graph_from_adjacency_matrix(traj.nn, mode = "undirected")
# is_connected(traj.graph) ## Ensure graph is connected
# 
# col.palette <- brewer.pal(length(levels(clusters)), "Set1")
# V(traj.graph)$size <- 1.0; V(traj.graph)$color <- sapply(as.integer(clusters), function(i) col.palette[[i]]);
# 
# l.use <- layout_with_kk(traj.graph)
# plot(traj.graph, layout = l.use, vertex.color = V(traj.graph)$color,
#      vertex.label = NA, edge.lty = "blank", vertex.frame.color = NA)

## Pick a k similar to the avg degree of our reference NN graph
k.use <- round(mean(degree(traj.graph)))
embeddings.knn <- lapply(embeddings, ComputeKNN, k = k.use)
# embeddings.snn <- lapply(embeddings, CalcSNN, k = k.use, prune.SNN = 0.0)

knn.simil <- sapply(embeddings.knn, function(knn.emb) {
  mean(sapply(1:ncol(knn.emb), function(i) CalcJaccard(knn.emb[,i], traj.nn[,i])))
})
print(knn.simil)

pdf("splatter_trajectory_neighborhood_score.pdf", width = 3.5, height = 4)
ggBarplot(knn.simil, fill.color = "skyblue")
dev.off()


## Calculate pairwise distances between trajectory sections
path.step.5 <- GetPathStep(metadata, step.size = 5, make.factor = T)
traj.dist <- CalcPairwiseDist(norm.counts[var.genes,], path.step.5, use.median = F)

embeddings.cor <- sapply(embeddings, function(emb) {
  emb.dist <- CalcPairwiseDist(emb, path.step.5, use.median = F)
  cor(traj.dist, emb.dist)
})
print(embeddings.cor)

pdf("splatter_trajectory_path_dist_cor.pdf", width = 3.5, height = 4)
ggBarplot(embeddings.cor, fill.color = "skyblue")
dev.off()

save.image("splatter_trajectory_analysis.RData")

#### Evaluate how SWNE visualization changes with k ####

setwd("/media/Scratch_SSD/Yan/R_Analysis/swne-dev/simulations")
load("splatter_trajectory_analysis.RData")
library(swne)

k.range <- c(3, 6, 9, 12)
for (k in k.range) {
  nmf.res <- RunNMF(norm.counts[var.genes,], k = k, alpha = 0, init = "random", 
                    n.cores = n.cores, loss = loss, n.rand.init = 5)
  swne.embedding <- EmbedSWNE(nmf.res$H, SNN = snn, alpha.exp = 2.5, snn.exp = 0.1, n_pull = 4,
                              pca.red = T, snn.factor.proj = T, dist.use = "cosine")
  swne.embedding$H.coords$name <- ""
  
  pdf(paste("splatter_trajectory_swne_k", k, ".pdf", sep = ""), width = 4, height = 4)
  print(PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = clusters, do.label = T, label.size = 4.5, 
                 pt.size = 1.5, seed = color.seed, show.legend = F))
  dev.off()
}
