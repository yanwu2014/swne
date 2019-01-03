## Read data
library(Seurat)
library(swne)
library(monocle)
source("simulations/splatter_helper_functions.R")

## Load monocle2 object
load("hemato_monocle_orig_cds.RData")

## Create Seurat object and create SNN
counts <- ReadData("hemato_expr_debatched.tsv", make.sparse = T)
info.genes <- scan("hemato_info_genes.txt", sep = "\n", what = character())

clusters.df <- read.table("hemato_cluster_mapping.csv", sep = ",")
clusters <- clusters.df[[2]]; names(clusters) <- clusters.df[[1]];

counts <- counts[,names(clusters)]
counts <- FilterData(counts, min.samples.frac = 0.005, trim = 0.005, min.nonzero.features = 200)
info.genes <- info.genes[info.genes %in% rownames(counts)]
dim(counts)

se.obj <- CreateSeuratObject(counts)
se.obj <- SetIdent(se.obj, cells.use = names(clusters), clusters)
se.obj@ident <- plyr::revalue(se.obj@ident,
                              c("1" = 'Ery', "2" = 'Ery', "3" = 'Ery', "4" = 'Ery', "5" = 'Ery',
                                "6" = 'Ery', "7" = 'MP/EP', "8" = 'MK', "9" = 'GMP', "10" = 'GMP',
                                "11" = 'DC', "12" = 'Bas', "13" = 'Bas', "14" = 'M', "15" = 'M',
                                "16" = 'Neu', "17" = 'Neu', "18" = 'Eos', "19" = 'lymphoid'))
rm(counts)

se.obj <- SubsetData(se.obj, ident.remove = c("lymphoid", "DC"))
se.obj@data <- ScaleCounts(se.obj@raw.data[,se.obj@cell.names], method = "log")
se.obj@scale.data <- se.obj@data - Matrix::rowMeans(se.obj@data)

se.obj <- RunPCA(se.obj, pc.genes = info.genes, do.print = F, pcs.compute = 40)
PCElbowPlot(se.obj, num.pc = 40)

se.obj <- RunTSNE(se.obj, dims.use = 1:15)
se.obj <- BuildSNN(se.obj, dims.use = 1:15, k.param = 40, prune.SNN = 0.0, force.recalc = T)

## Run UMAP
se.obj <- RunUMAP(se.obj, reduction.use = "pca", dims.use = 1:15)

## Pull out clusters and Set cluster colors
clusters <- se.obj@ident; names(clusters) <- se.obj@cell.names;
cluster_colors <- c("Bas" = "#ff6347", "Eos" = "#EFAD1E",
                    "Ery" = "#8CB3DF", "M" = "#53C0AD", "MP/EP" = "#4EB859",
                    "GMP" = "#D097C4", "MK" = "#ACC436", "Neu" = "#F5918A")

## SWNE parameters
k.range <- seq(2,20,2)
n.cores <- 24

norm.counts <- se.obj@data[,names(clusters)]
k.err <- FindNumFactors(norm.counts[info.genes,], k.range = k.range, n.cores = n.cores,
                        seed = 32590, do.plot = F)

library(ggplot2)
pdf("hemato_mdl_selection.pdf", width = 5, height = 4)
PlotFactorSelection(k.err, font.size = 15) +
  theme(axis.title = element_blank(), legend.title = element_blank())
dev.off()


k <- 12
nmf.res <- RunNMF(norm.counts[info.genes,], k = k, init = "ica", n.cores = n.cores, ica.fast = F)
nmf.res$W <- ProjectFeatures(norm.counts, nmf.res$H, n.cores = n.cores)
nmf.scores <- nmf.res$H

swne.embedding <- EmbedSWNE(nmf.scores, SNN = se.obj@snn, alpha.exp = 1.5, snn.exp = 1, n_pull = 3,
                            proj.method = "sammon", dist.use = "cosine")

swne.embedding <- RenameFactors(swne.embedding, name.mapping = c("factor_4" = "Erythrocyte differentiation",
                                                                 "factor_5" = "Metal binding",
                                                                 "factor_9" = "Epigenetic regulation",
                                                                 "factor_10" = "HSC maintenance",
                                                                 "factor_11" = "Inflammation"))
## Embed selected genes onto swne plot
genes.embed <- c("Apoe", "Mt2", "Gpr56", "Sun2", "Flt3")
swne.embedding <- EmbedFeatures(swne.embedding, nmf.res$W, genes.embed, n_pull = 3)

## SWNE plots
pdf("hemato_swne_plot.pdf", width = 5.5, height = 5)
PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = clusters, do.label = T,
         label.size = 3.5, pt.size = 2, show.legend = F) +
  scale_color_manual(values = cluster_colors)
dev.off()

pdf("hemato_swne_plot_nolabels.pdf", width = 5.5, height = 5)
PlotSWNE(swne.embedding, alpha.plot = 0.5, sample.groups = clusters, do.label = F,
         label.size = 0, pt.size = 2, show.legend = F) +
  scale_color_manual(values = cluster_colors)
dev.off()

## tSNE plots
tsne.scores <- GetCellEmbeddings(se.obj, reduction.type = "tsne")
pdf("hemato_tsne_plot.pdf", width = 5, height = 5)
PlotDims(tsne.scores, sample.groups = clusters, show.legend = F, show.axes = F,
         alpha.plot = 0.5, label.size = 3.5, pt.size = 1.5) +
  scale_color_manual(values = cluster_colors)
dev.off()

pdf("hemato_tsne_plot_nolabels.pdf", width = 5, height = 5)
PlotDims(tsne.scores, sample.groups = clusters, show.legend = F, show.axes = F,
         alpha.plot = 0.5, pt.size = 1.5, label.size = 0) +
  scale_color_manual(values = cluster_colors)
dev.off()


## UMAP plots
umap.scores <- GetCellEmbeddings(se.obj, reduction.type = "umap")
pdf("hemato_umap_plot.pdf", width = 5, height = 5)
PlotDims(umap.scores, sample.groups = clusters, show.legend = F, show.axes = F,
         alpha.plot = 0.5, label.size = 3.5, pt.size = 1.5) +
  scale_color_manual(values = cluster_colors)
dev.off()

pdf("hemato_umap_plot_nolabels.pdf", width = 5, height = 5)
PlotDims(umap.scores, sample.groups = clusters, show.legend = F, show.axes = F,
         alpha.plot = 0.5, pt.size = 1.5, label.size = 0) +
  scale_color_manual(values = cluster_colors)
dev.off()



## Overlay SWNE plot with pseudotime
pseudotime <- cds$Pseudotime; names(pseudotime) <- colnames(cds@reducedDimS);

pdf("hemato_swne_pseudotime_plot.pdf", width = 6, height = 5)
FeaturePlotSWNE(swne.embedding, pseudotime, alpha.plot = 0.4, label.size = 3.5, pt.size = 1.5)
dev.off()

pdf("hemato_swne_pseudotime_plot_nolabels.pdf", width = 6, height = 5)
FeaturePlotSWNE(swne.embedding, pseudotime, alpha.plot = 0.4, label.size = 0, pt.size = 1.5)
dev.off()

pdf("hemato_tsne_pseudotime_plot.pdf", width = 5.5, height = 5)
FeaturePlotDims(tsne.scores, pseudotime, alpha.plot = 0.4, show.axes = F, pt.size = 1.5)
dev.off()

pdf("hemato_umap_pseudotime_plot.pdf", width = 5.5, height = 5)
FeaturePlotDims(umap.scores, pseudotime, alpha.plot = 0.4, show.axes = F, pt.size = 1.5)
dev.off()


## Identify factors using the gene loadings matrix
factor.genes.df <- SummarizeAssocFeatures(nmf.res$W, features.return = 8)
write.table(factor.genes.df, file = "hemato_factor_markers.tsv", sep = "\t")

## Identify cluster marker genes
cluster.genes.df <- FindAllMarkers(se.obj, logfc.threshold = 0.25, only.pos = T)
top.cluster.genes.df <- Reduce("rbind", by(cluster.genes.df, cluster.genes.df$cluster, head, n = 2))

cluster.genes.heat <- t(apply(norm.counts[genes.embed,], 1, function(x) {
  x <- (x - mean(x)); tapply(x, clusters, mean);
}))

pdf("hemato_marker_gene_heatmap.pdf", width = 5, height = 3)
ggHeat(cluster.genes.heat, clustering = "both", x.lab.size = 16, y.lab.size = 15)
dev.off()

save.image("hemato_swne.RData")



## Quantitative evaluation of embeddings

## Compute Embeddings
library(destiny)
dm <- DiffusionMap(t(as.matrix(norm.counts[info.genes,])), k = 20, n_eigs = 2)
diff.scores <- dm@eigenvectors; rownames(diff.scores) <- colnames(norm.counts);

library(RDRToolbox)
isomap.scores <- Isomap(t(as.matrix(norm.counts[info.genes,])), dims = 2, k = 20)$dim2
rownames(isomap.scores) <- colnames(norm.counts)

lle.scores <- LLE(t(as.matrix(norm.counts[info.genes,])), dim = 2, k = 20)
rownames(lle.scores) <- colnames(norm.counts)

mds.scores <- cmdscale(dist(t(as.matrix(norm.counts[info.genes,]))), k = 2)
pc.scores <- GetCellEmbeddings(se.obj, dims.use = 1:2)

embeddings <- list(swne = t(as.matrix(swne.embedding$sample.coords)),
                   tsne = t(tsne.scores),
                   pca = t(pc.scores),
                   lle = t(lle.scores),
                   mds = t(mds.scores),
                   isomap = t(isomap.scores),
                   dmaps = t(diff.scores),
                   umap = t(umap.scores))

## Calculate neighborhood fidelity
n.neighbors <- 40
ref.knn <- ComputeKNN(norm.counts[info.genes,], k = n.neighbors)

## Compute kNN for embeddings
embeddings.knn <- lapply(embeddings, ComputeKNN, k = n.neighbors)

knn.simil <- sapply(embeddings.knn, function(knn.emb) {
  mean(sapply(1:ncol(knn.emb), function(i) CalcJaccard(knn.emb[,i], ref.knn[,i])))
})
print(knn.simil)

## Compute cluster-time distance correlations
metadata.df <- data.frame(Group = clusters, Step = order(pseudotime[names(clusters)]))
clusters.steps <- GetPathStep(metadata.df, step.size = 50, make.factor = T)
traj.dist <- CalcPairwiseDist(norm.counts[info.genes,], clusters.steps)

embeddings.cor <- sapply(embeddings, function(emb) {
  emb.dist <- CalcPairwiseDist(emb, clusters.steps)
  cor(traj.dist, emb.dist)
})
print(embeddings.cor)


library(ggplot2)
library(ggrepel)

pdf("hemato_quant_eval.pdf", width = 6.5, height = 4.5)
scatter.df <- data.frame(x = knn.simil, y = embeddings.cor, name = names(embeddings))
ggplot(scatter.df, aes(x, y)) + geom_point(size = 2, alpha = 1) +
  theme_classic() + theme(legend.position = "none", text = element_text(size = 16)) +
  xlab("Neighborhood Score") + ylab("Path-Time Distance Correlation") +
  geom_text_repel(aes(x, y, label = name), size = 5) +
  xlim(0, max(knn.simil)) + ylim(0, max(embeddings.cor))
dev.off()

save.image("hemato_swne.RData")


## Create original monocle2 plots
rge.scores <- t(cds@reducedDimS)[names(clusters), 1:2]

pdf("hemato_rge_2d_plot.pdf", width = 5, height = 5)
PlotDims(rge.scores, sample.groups = clusters, show.legend = F, show.axes = F,
         alpha.plot = 0.5, label.size = 7, pt.size = 1.5) +
  scale_color_manual(values = cluster_colors)
dev.off()

pdf("hemato_rge_2d_pseudotime_plot.pdf", width = 5.5, height = 5)
FeaturePlotDims(rge.scores, pseudotime, alpha.plot = 0.4, show.axes = F, pt.size = 1.5)
dev.off()

cds <- reduceDimension(cds, max_components = 10, norm_method = "log", reduction_method = "DDRTree")
cds <- cds[,names(clusters)]
cds$cell_type <- clusters

pdf("hemato_rge_tree_plot.pdf", width = 5, height = 5)
plot_complex_cell_trajectory(cds, color_by = "cell_type") +
  scale_color_manual(values = cluster_colors) +
  theme(legend.text = element_text(size = 16), legend.title = element_blank())
dev.off()

cds$Pseudotime <- pseudotime[names(clusters)]
pdf("hemato_rge_tree_pseudotime_plot.pdf", width = 5, height = 5)
plot_complex_cell_trajectory(cds, color_by = "Pseudotime") +
  scale_color_distiller(palette = "YlOrRd", direction = 1,
                        guide_colorbar(title = "Pseudotime", label = T)) +
  theme(legend.text = element_text(size = 12, angle = 90, vjust = 0.7))
dev.off()

## Diffusion maps plots
pdf("hemato_dmaps_plot.pdf", width = 5, height = 5)
PlotDims(diff.scores, sample.groups = clusters, show.legend = F, show.axes = F,
         alpha.plot = 0.5, label.size = 7, pt.size = 1.5) +
  scale_color_manual(values = cluster_colors)
dev.off()

pdf("hemato_dmaps_pseudotime_plot.pdf", width = 5.5, height = 5)
FeaturePlotDims(diff.scores, pseudotime, alpha.plot = 0.4, show.axes = F, pt.size = 1.5)
dev.off()


## Validate gene embeddings
load("hemato_swne.RData")
library(swne)

gene <- "Flt3"
gene.swne.embedding <- swne.embedding
gene.swne.embedding <- EmbedFeatures(gene.swne.embedding, nmf.res$W, unique(c(genes.embed, gene)),
                                     n_pull = 6, scale.cols = T)
gene.swne.embedding$feature.coords <- subset(gene.swne.embedding$feature.coords, name == gene)
gene.swne.embedding$H.coords$name <- ""

pdf(paste0("hemato_", gene, "_feature_plot.pdf"), width = 4.5, height = 4)
FeaturePlotSWNE(gene.swne.embedding, norm.counts[gene,], pt.size = 1, label.size = 0)
dev.off()


pdf("hemato_gene_embedding_check.pdf", width = 6, height = 5)
gene.logfc.df <- CheckGeneEmbedding(nmf.res$W, norm.counts, genes.embed, clusters, min.cluster.logfc = 1.5,
                                    min.factor.logfc = 2.5, label.size = 6.5, font.size = 18)
dev.off()

## Find some bad genes to embed to see what happens
gene.logfc.df <- gene.logfc.df[order(rowMeans(gene.logfc.df[,1:2])),]
head(gene.logfc.df)

save.image("hemato_swne.RData")


## Systematically see where diff and non-diff genes tend to be embedded
all.diff.genes <- rownames(subset(gene.logfc.df, cluster > 1.5 & factor > 2.5))
non.diff.genes <- rownames(subset(gene.logfc.df, cluster < 1.5 | factor < 2.5))

diff.gene.swne.embedding <- swne.embedding
diff.gene.swne.embedding <- EmbedFeatures(diff.gene.swne.embedding, nmf.res$W, all.diff.genes, n_pull = 6)
# diff.gene.swne.embedding <- EmbedFeatures(diff.gene.swne.embedding, nmf.res$W, non.diff.genes, n_pull = 6)

feature.emb <- diff.gene.swne.embedding$feature.coords
diff.gene.swne.embedding$feature.coords <- NULL

sample.coords.df <- diff.gene.swne.embedding$sample.coords
sample.coords.df$sample.groups <- clusters

H.coords.df <- subset(diff.gene.swne.embedding$H.coords, name != "")

pdf("hemato_diff_gene_embeddings.pdf", width = 4.5, height = 4)
ggplot() +
  geom_point(data = sample.coords.df, aes(x, y),
             alpha = 0.5, size = 1.25) +
  theme_classic() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                          axis.ticks = element_blank(), axis.line = element_blank(),
                          axis.text = element_blank(), legend.text = element_text(size = 12),
                          legend.title = element_blank()) +
  # geom_point(data = H.coords.df, aes(x, y), size = 2.5, color = "darkblue") +
  geom_hex(data = feature.emb, mapping = aes(x, y), bins = 20, alpha = 0.5) +
  scale_fill_gradient(low = "skyblue", high = "darkred")
dev.off()


## Make some histograms of the data
pdf("hemato_matrix_hist.pdf", width = 3, height = 3)
hist(as.matrix(norm.counts), xlab = "Scaled expression", main = NULL)
dev.off()

pdf("hemato_nonzero_matrix_hist.pdf", width = 3, height = 3)
hist(norm.counts@x, xlab = "Scaled expression", main = NULL)
dev.off()
