library(Seurat)
library(swne)
source("simulations/splatter_helper_functions.R")

## Load data
counts <- ReadData("snDropSeq_vc_cer_counts.tsv")
metadata.df <- read.table("snDropSeq_metadata.txt", sep = "\t", header = T, stringsAsFactors = F)
rownames(metadata.df) <- paste(metadata.df$Identity, rownames(metadata.df), sep = "_")

## Filter counts
counts <- FilterData(counts, min.samples.frac = 0.001, trim = 0.001, min.nonzero.features = 200)
metadata.df <- metadata.df[colnames(counts),]
dim(counts)

## Create Seurat object and normalize data
se.obj <- CreateSeuratObject(counts,  meta.data = metadata.df)
se.obj@data <- ScaleCounts(se.obj@raw.data[,se.obj@cell.names], method = "log", adj.var = T)
rm(counts)

## Select overdispersed genes
n.var.genes <- 3e3
var.df <- AdjustVariance(se.obj@raw.data, verbose = F)
var.genes <- rownames(var.df[order(var.df$lp),])[1:n.var.genes]
se.obj@var.genes <- var.genes

## Center data and run PCA
se.obj@scale.data <- se.obj@data[var.genes,] - Matrix::rowMeans(se.obj@data[var.genes,])
se.obj <- RunPCA(se.obj, pcs.compute = 40, do.print = F)
PCElbowPlot(se.obj, num.pc = 40)

pcs.use <- 30
se.obj <- BuildSNN(se.obj, dims.use = 1:pcs.use, k = 30, prune.SNN = 0.05, force.recalc = T)
se.obj <- RunTSNE(se.obj, dims.use = 1:pcs.use)
se.obj <- RunUMAP(se.obj, dims.use = 1:pcs.use)

clusters <- se.obj@ident; names(clusters) <- se.obj@cell.names; levels(clusters);
clusters <- plyr::revalue(clusters, replace = 
                            c("Ex1" = "Ex_L2/3", "Ex3a" = "Ex_L4", "Ex3b" = "Ex_L4", "Ex3c" = "Ex_L4", "Ex3d" = "Ex_L4",
                              "Ex4" = "Ex_L4/5", "Ex5a" = "Ex_L5", "Ex5b" = "Ex_L5", "Ex6a" = "Ex_L6", "Ex6b" = "Ex_L6", 
                              "Ex8" = "Ex_L6b", "In1a" = "In1", "In1b" = "In1", "In1c" = "In1", "In4a" = "In4",
                              "In4b" = "In4", "In6a" = "In6", "In6b" = "In6", "In7" = "In7/8",
                              "In8" = "In7/8"))
levels(clusters)

## Region of origin
region <- factor(metadata.df$Brain_Region); names(region) <- se.obj@cell.names;
levels(region) <- c("Visual Cortex", "Lateral Cerebellum")

## Normalize counts matrix
norm.counts <- se.obj@data

#### Run SWNE on full dataset ####
k.range <- seq(2,30,2)
n.cores <- 24
seed <- 223464

## Unguided NMF
k.err <- FindNumFactors(norm.counts[var.genes,], k.range = k.range, n.cores = n.cores, do.plot = T,
                        seed = seed)

library(ggplot2)
pdf("snDropSeq_mdl_selection.pdf", width = 5, height = 4)
PlotFactorSelection(k.err, font.size = 15) + 
  theme(legend.title = element_blank(), axis.title = element_blank())
dev.off()

k <- 28
nmf.res <- RunNMF(norm.counts[var.genes,], k = k, init = "ica", n.cores = n.cores, ica.fast = T)
nmf.res$W <- ProjectFeatures(norm.counts, nmf.res$H, n.cores = n.cores)
nmf.scores <- nmf.res$H

swne.embedding <- EmbedSWNE(nmf.scores, se.obj@snn, alpha.exp = 1.5, snn.exp = 0.1, 
                            n_pull = 3, dist.use = "cosine")

## Rename NMFs
swne.embedding <- RenameFactors(swne.embedding, name.mapping =
                                  c("factor_11" = "Axon projection", "factor_12" = "Glutamate transport",
                                    "factor_21" = "Myelin formation","factor_26" = "Cell junctions"))

## Embed genes
genes.embed <- c("PLP1", "CBLN2", "LHFPL3", "SLC1A2", "FSTL5", "NRGN", "GRIK1", 
                 "NTNG1", "DAB1", "HS3ST2", "DCC", "POSTN")
swne.embedding <- EmbedFeatures(swne.embedding, nmf.res$W, genes.embed, n_pull = 3, scale.cols = F)

plot.seed <- 23529080
pdf("snDropSeq_swne_plot_clusters.pdf", width = 6.5, height = 6.5)
PlotSWNE(swne.embedding, alpha.plot = 0.3, sample.groups = clusters, do.label = T, 
         label.size = 3.5, pt.size = 0.75, show.legend = F, seed = plot.seed)
dev.off()

pdf("snDropSeq_swne_plot_nolabel.pdf", width = 6.5, height = 6.5)
PlotSWNE(swne.embedding, alpha.plot = 0.3, sample.groups = clusters, do.label = T, 
         label.size = 0, pt.size = 0.75, show.legend = F, seed = plot.seed)
dev.off()

tsne.scores <- GetCellEmbeddings(se.obj, reduction.type = "tsne")
pdf("snDropSeq_tsne_plot_clusters.pdf", width = 6.5, height = 6.5)
PlotDims(tsne.scores, sample.groups = clusters, pt.size = 0.5, label.size = 4, 
         alpha.plot = 0.4, show.legend = F, show.axes = F, seed = plot.seed)
dev.off()

pdf("snDropSeq_tsne_plot_nolabel.pdf", width = 6.5, height = 6.5)
PlotDims(tsne.scores, sample.groups = clusters, pt.size = 0.5, label.size = 0, 
         alpha.plot = 0.4, show.legend = F, show.axes = F, seed = plot.seed)
dev.off()


## UMAP plots
umap.scores <- GetCellEmbeddings(se.obj, reduction.type = "umap")
pdf("snDropSeq_umap_plot.pdf", width = 6.5, height = 6.5)
PlotDims(umap.scores, sample.groups = clusters, show.legend = F, show.axes = F, 
         alpha.plot = 0.4, label.size = 4, pt.size = 0.5, seed = plot.seed)
dev.off()

pdf("snDropSeq_umap_plot_nolabels.pdf", width = 6.5, height = 6.5)
PlotDims(umap.scores, sample.groups = clusters, show.legend = F, show.axes = F, 
         alpha.plot = 0.4, pt.size = 0.5, label.size = 0, seed = plot.seed)
dev.off()


## Associate factors with genes using the gene loadings (W) matrix
factor.genes.df <- SummarizeAssocFeatures(nmf.res$W, features.return = 8)
write.table(factor.genes.df, file = "snDropSeq_factor_markers.tsv", sep = "\t")

## Find cluster markers
se.obj <- SetIdent(se.obj, cells.use = names(clusters), ident.use = clusters)
cluster.genes.df <- FindAllMarkers(se.obj, logfc.threshold = 0.5, only.pos = T)

top.cluster.genes.df <- Reduce("rbind", by(cluster.genes.df, cluster.genes.df$cluster, head, n = 1))
# cluster.genes.plot <- unique(c(top.cluster.genes.df$gene, genes.embed))
cluster.genes.plot <- genes.embed

cluster.genes.heat <- norm.counts[cluster.genes.plot,]
cluster.genes.heat <- cluster.genes.heat - Matrix::rowMeans(cluster.genes.heat)
cluster.genes.heat <- t(apply(cluster.genes.heat, 1, function(x) tapply(x, clusters, mean)))

pdf("snDropSeq_gene_cluster_heatmap.pdf", width = 7.5, height = 5)
ggHeat(cluster.genes.heat, clustering = "both", x.lab.size = 16, y.lab.size = 15)
dev.off()


## Quantitative evaluations
pc.scores <- GetCellEmbeddings(se.obj, dims.use = 1:2)
clust.dist <- CalcPairwiseDist(norm.counts[var.genes,], clusters)

embeddings <- list(swne = t(as.matrix(swne.embedding$sample.coords)),
                   tsne = t(tsne.scores),
                   pca = t(pc.scores),
                   # lle = t(lle.scores),
                   # mds = t(mds.scores),
                   # isomap = t(isomap.scores),
                   # dmaps = t(diff.scores),
                   umap = t(umap.scores))

## Calculate neighborhood fidelity
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

pdf("snDropSeq_quant_eval.pdf", width = 6, height = 5)
ggplot(scatter.df, aes(x, y)) + geom_point(size = 2, alpha = 1) +
  theme_classic() + theme(legend.position = "none", text = element_text(size = 16)) +
  xlab("Silhouette Score") + ylab("Cluster Distance Correlation") + 
  geom_text_repel(aes(x, y, label = name), size = 5) + 
  xlim(0, max(si.avg.scores)) + ylim(0, max(embeddings.cor))
dev.off()


## Validate gene embeddings
load("snDropSeq_swne.RData")
library(swne)

gene <- "CADM2"
gene.swne.embedding <- swne.embedding
gene.swne.embedding <- EmbedFeatures(gene.swne.embedding, nmf.res$W, unique(c(genes.embed, gene)), 
                                     n_pull = 3, scale.cols = T)
gene.swne.embedding$feature.coords <- subset(gene.swne.embedding$feature.coords, name == gene)
gene.swne.embedding$H.coords$name <- ""

pdf(paste0("snDropSeq_", gene, "_feature_plot_BuPu.pdf"), width = 4.5, height = 4)
FeaturePlotSWNE(gene.swne.embedding, norm.counts[gene,], pt.size = 0.35, alpha= 0.3, 
                label.size = 0, color.palette = "BuPu") + 
  theme(legend.text = element_blank())
dev.off()


pdf("snDropSeq_gene_embedding_check.pdf", width = 6, height = 5)
gene.logfc.df <- CheckGeneEmbedding(nmf.res$W, norm.counts, genes.embed, clusters, min.cluster.logfc = 1.5,
                                    min.factor.logfc = 2, label.size = 6, font.size = 18)
dev.off()

## Find some bad genes to embed to see what happens
gene.logfc.df <- gene.logfc.df[order(rowMeans(gene.logfc.df[,1:2])),]
head(gene.logfc.df)

save.image("snDropSeq_swne.RData")


## Systematically see where diff and non-diff genes tend to be embedded
all.diff.genes <- rownames(subset(gene.logfc.df, cluster > 1.5 & factor > 2))
non.diff.genes <- rownames(subset(gene.logfc.df, cluster < 1.5 | factor < 2))

diff.gene.swne.embedding <- swne.embedding
diff.gene.swne.embedding <- EmbedFeatures(diff.gene.swne.embedding, nmf.res$W, all.diff.genes, n_pull = 3)
# diff.gene.swne.embedding <- EmbedFeatures(diff.gene.swne.embedding, nmf.res$W, non.diff.genes, n_pull = 3)

feature.emb <- diff.gene.swne.embedding$feature.coords
diff.gene.swne.embedding$feature.coords <- NULL

sample.coords.df <- diff.gene.swne.embedding$sample.coords
sample.coords.df$sample.groups <- clusters
H.coords.df <- subset(diff.gene.swne.embedding$H.coords, name != "")

pdf("snDropSeq_diff_gene_embeddings_BuPu.pdf", width = 4.5, height = 4)
ggplot() + 
  geom_point(data = sample.coords.df, aes(x, y),
             alpha = 1, size = 1.25, color = "grey") +
  theme_classic() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                          axis.ticks = element_blank(), axis.line = element_blank(),
                          axis.text = element_blank(), legend.text = element_blank(),
                          legend.title = element_blank()) + 
  geom_hex(data = feature.emb, mapping = aes(x, y), bins = 25, alpha = 0.65) + 
  scale_fill_gradient(low = "skyblue", high = "darkviolet")
dev.off()



#### Project C1 data onto snDropSeq data ####
c1.obj <- readRDS("C1_cortex_seurat.Robj")
c1.norm.counts <- ScaleCounts(as(as.matrix(c1.obj@raw.data[,c1.obj@cell.names]), "dgCMatrix"), method = "ft")

c1.clusters <- c1.obj@ident; names(c1.clusters) <- c1.obj@cell.names; levels(c1.clusters);
c1.clusters <- plyr::revalue(c1.clusters, replace = 
                               c("Ex1" = "Ex_L2/3", "Ex2" = "Ex_L4", "Ex3" = "Ex_L4", 
                                 "Ex4" = "Ex_L4/5", "Ex5" = "Ex_L5", "Ex6" = "Ex_L6", 
                                 "Ex7" = "Ex_L6", "Ex8" = "Ex_L6b", 
                                 "In7" = "In7/8", "In8" = "In7/8"))
levels(c1.clusters)

genes.project <- intersect(var.genes, rownames(c1.norm.counts))
c1.nmf.scores <- ProjectSamples(c1.norm.counts, nmf.res$W, features.use = genes.project, n.cores = n.cores)

c1.snn <- ProjectSNN(c1.norm.counts, norm.counts, n.pcs = 30, features.use = genes.project, k = 40)
c1.swne.embedding <- ProjectSWNE(swne.embedding, c1.nmf.scores, SNN = c1.snn, 
                                 alpha.exp = 1.5, snn.exp = 0.1, n_pull = 3)

merged.embedding <- swne.embedding
merged.embedding$sample.coords <- rbind(merged.embedding$sample.coords, c1.swne.embedding$sample.coords)
merged.clusters <- c(as.character(clusters), as.character(c1.clusters))
names(merged.clusters) <- c(names(clusters), names(c1.clusters))

neurons <- names(merged.clusters[grepl("Ex|In", merged.clusters)])
merged.clusters <- factor(merged.clusters)

cluster.colors <- ExtractSWNEColors(swne.embedding, clusters, seed = plot.seed)
cluster.colors <- cluster.colors[grepl("Ex|In", names(cluster.colors))]
cluster.colors[["In5"]] <- "#00BFFF"

pdf("snDropSeq_swne_plot_c1_neurons.pdf", width = 6.5, height = 6.5)
PlotSWNE(merged.embedding, alpha.plot = 0.5, sample.groups = merged.clusters, do.label = T, 
         pt.size = 1.5, samples.plot = names(c1.clusters), show.legend = F, seed = plot.seed) +
  scale_color_manual(values = cluster.colors)
dev.off()

pdf("snDropSeq_swne_plot_c1_neurons_nolabels.pdf", width = 6.5, height = 6.5)
PlotSWNE(merged.embedding, alpha.plot = 0.5, sample.groups = merged.clusters, do.label = T, 
         pt.size = 1.5, samples.plot = names(c1.clusters), show.legend = F, 
         seed = plot.seed, label.size = 0) +
  scale_color_manual(values = cluster.colors)
dev.off()


tech <- as.character(merged.clusters); names(tech) <- names(merged.clusters);
tech[names(clusters)] <- "snDropSeq"
tech[names(c1.clusters)] <- "C1"
tech <- factor(tech, levels = c("C1", "snDropSeq"))

merged.embedding$sample.coords$tech <- tech
merged.embedding$sample.coords <- merged.embedding$sample.coords[order(merged.embedding$sample.coords$tech, decreasing = T),]
tech <- tech[rownames(merged.embedding$sample.coords)]

pdf("snDropSeq_swne_plot_neurons_tech.pdf", width = 6.5, height = 6.5)
PlotSWNE(merged.embedding, alpha.plot = 0.5, sample.groups = tech, do.label = F, 
         pt.size = 0.8, samples.plot = rownames(merged.embedding$sample.coords), 
         show.legend = F, seed = plot.seed, label.size = 0)
dev.off()

save.image("snDropSeq_swne.RData")


## Make some histograms of the data
pdf("snDropSeq_matrix_hist.pdf", width = 3, height = 3)
hist(as.matrix(norm.counts), xlab = "Scaled expression", main = NULL)
dev.off()

pdf("snDropSeq_nonzero_matrix_hist.pdf", width = 3, height = 3)
hist(norm.counts@x, xlab = "Scaled expression", main = NULL)
dev.off()