library(Seurat)
library(swne)

## Load data
counts <- ReadData("snDropSeq_vc_cer_counts.tsv")
metadata.df <- read.table("snDropSeq_metadata.txt", sep = "\t", header = T, stringsAsFactors = F)
rownames(metadata.df) <- paste(metadata.df$Identity, rownames(metadata.df), sep = "_")

counts <- FilterData(counts, min.samples.frac = 0.001, trim = 0.001, min.nonzero.features = 200)
metadata.df <- metadata.df[colnames(counts),]
dim(counts)

## Use Seurat for PCA, tSNE, SNN
se.obj <- CreateSeuratObject(counts, normalization.method = "LogNormalize", scale.factor = median(Matrix::colSums(counts)),
                             meta.data = metadata.df)
se.obj <- FindVariableGenes(se.obj, x.low.cutoff = 0.01, y.cutoff = 0.5)
length(se.obj@var.genes)
rm(counts)

se.obj <- ScaleData(se.obj, genes.use = se.obj@var.genes, model.use = "negbinom", 
                    vars.to.regress = c("nUMI"))
se.obj <- RunPCA(se.obj, pcs.compute = 40, do.print = F)
PCElbowPlot(se.obj, num.pc = 40)

pcs.use <- 20
se.obj <- RunTSNE(se.obj, dims.use = 1:pcs.use)

clusters <- se.obj@ident; names(clusters) <- se.obj@cell.names; levels(clusters);
clusters <- plyr::revalue(clusters, replace = 
                            c("Ex1" = "Ex", "Ex3a" = "Ex", "Ex3b" = "Ex", "Ex3c" = "Ex", "Ex3d" = "Ex",
                              "Ex4" = "Ex", "Ex5a" = "Ex", "Ex5b" = "Ex", "Ex6a" = "Ex", "Ex6b" = "Ex",
                              "Ex8" = "Ex", "In1a" = "In1", "In1b" = "In1", "In1c" = "In1", "In4a" = "In4",
                              "In4b" = "In4", "In6a" = "In6", "In6b" = "In6", "In7" = "In7/8",
                              "In8" = "In7/8"))
levels(clusters)

## Region of origin
region <- factor(metadata.df$Brain_Region); names(region) <- se.obj@cell.names;
levels(region) <- c("Visual Cortex", "Lateral Cerebellum")

## Normalize counts matrix
norm.counts <- ScaleCounts(se.obj@raw.data[,se.obj@cell.names], batch = NULL, method = "ft", adj.var = T)

## NMF analysis
loss <- "mse"
k.range <- seq(2,24,2)
n.cores <- 24
seed <- 223464

## Unguided NMF
k.res <- FindNumFactors(norm.counts[se.obj@var.genes,], k.range = k.range, n.cores = n.cores, do.plot = T,
                        na.frac = 0.25, seed = 2530980, loss = loss)
k.res$k

k <- 20
nmf.res <- RunNMF(norm.counts[se.obj@var.genes,], k = k, alpha = 0, init = "ica", n.cores = n.cores,
                  loss = loss, init.zeros = "random")
nmf.res$W <- ProjectFeatures(norm.counts, nmf.res$H, loss = loss, n.cores = n.cores)
nmf.scores <- nmf.res$H

se.obj <- BuildSNN(se.obj, dims.use = 1:pcs.use, k = 30, prune.SNN = 1/15, force.recalc = T)
swne.embedding <- EmbedSWNE(nmf.scores, se.obj@snn, alpha.exp = 1.75, snn.exp = 0.1, 
                            n_pull = 6, dist.use = "cosine")

## Rename NMFs
swne.embedding <- RenameFactors(swne.embedding, name.mapping =
                                  c("factor_2" = "Myelin", "factor_9" = "Immune response", 
                                    "factor_18" = "Cell junctions"))

## Embed genes
genes.embed <- c("PLP1", "CBLN2", "LHFPL3", "SLC1A2", "FSTL5", "NRGN", "GRIK1")
swne.embedding <- EmbedFeatures(swne.embedding, nmf.res$W, genes.embed, n_pull = 4, scale.cols = F)

pdf("snDropSeq_swne_plot_clusters.pdf", width = 6.5, height = 6.5)
PlotSWNE(swne.embedding, alpha.plot = 0.3, sample.groups = clusters, do.label = T, 
         label.size = 3.5, pt.size = 1.0, show.legend = F, seed = seed)
dev.off()

pdf("snDropSeq_swne_plot_nolabel.pdf", width = 6.5, height = 6.5)
PlotSWNE(swne.embedding, alpha.plot = 0.3, sample.groups = clusters, do.label = T, 
         label.size = 0, pt.size = 1.0, show.legend = F, seed = seed)
dev.off()


tsne.scores <- GetCellEmbeddings(se.obj, reduction.type = "tsne")
pdf("snDropSeq_tsne_plot_clusters.pdf", width = 6, height = 6)
PlotDims(tsne.scores, sample.groups = clusters, pt.size = 0.5, label.size = 4, 
         alpha.plot = 0.4, show.legend = F, show.axes = F, seed = seed)
dev.off()

pdf("snDropSeq_tsne_plot_nolabel.pdf", width = 6, height = 6)
PlotDims(tsne.scores, sample.groups = clusters, pt.size = 0.5, label.size = 0, 
         alpha.plot = 0.4, show.legend = F, show.axes = F, seed = seed)
dev.off()


## Associate factors with cell clusters
clusters.list <- UnflattenGroups(clusters)
clusters.matrix <- t(swne:::.genesets_indicator(clusters.list, inv = F, return.numeric = T))
cluster.nmf.assoc <- FactorAssociation(clusters.matrix, nmf.scores, n.cores = n.cores, metric = "IC")

pdf("snDropSeq_cluster_factor_heatmap.pdf", width = 7.5, height = 4.5)
ggHeat(cluster.nmf.assoc, clustering = "both", x.lab.size = 14, y.lab.size = 12)
dev.off()

## Associate factors with genes using the gene loadings (W) matrix
gene.loadings <- nmf.res$W
top.factor.genes.df <- SummarizeAssocFeatures(gene.loadings, features.return = 6)
# write.table(top.factor.genes.df, file = "snDropSeq_factor_markers.tsv", sep = "\t")

## Make gene factor association heatmaps
gene.factor.heat <- gene.loadings[unique(top.factor.genes.df$feature),]

pdf("snDropSeq_gene_factor_heatmap.pdf", width = 7.5, height = 7.0)
ggHeat(gene.factor.heat, clustering = "both", x.lab.size = 14, y.lab.size = 12)
dev.off()

## Find cluster markers
se.obj <- SetIdent(se.obj, cells.use = names(clusters), ident.use = clusters)
cluster.genes.df <- FindAllMarkers(se.obj, logfc.threshold = 0.2, only.pos = T)
cluster.genes.mat <- UnflattenDataframe(cluster.genes.df, output.name = "avg_logFC",
                                        row.col = "gene", col.col = "cluster")
cluster.genes.mat[is.na(cluster.genes.mat)] <- 0

top.cluster.genes.df <- Reduce("rbind", by(cluster.genes.df, cluster.genes.df$cluster, head, n = 2))
cluster.genes.heat <- cluster.genes.mat[unique(top.cluster.genes.df$gene),]

pdf("snDropSeq_gene_cluster_heatmap.pdf", width = 6.5, height = 6)
ggHeat(cluster.genes.heat, clustering = "both", x.lab.size = 14, y.lab.size = 12)
dev.off()


## Validate gene embeddings
gene <- "CBLN2"
gene.swne.embedding <- swne.embedding
gene.swne.embedding$H.coords$name <- ""
gene.swne.embedding$feature.coords <- subset(gene.swne.embedding$feature.coords, name == gene)

pdf("snDropSeq_CBLN2_feature_plot.pdf", width = 4.5, height = 4.5)
FeaturePlotSWNE(gene.swne.embedding, norm.counts[gene,], pt.size = 0.5, label.size = 0)
dev.off()

save.image("snDropSeq_swne.RData")

#### Subset to Excitatory neurons to look at layer specificity ####
ex.se.obj <- SubsetData(se.obj, ident.use = c("Ex1", "Ex3a", "Ex3b", "Ex3c", "Ex3d", "Ex4", 
                                              "Ex5a", "Ex5b", "Ex6a", "Ex6b", "Ex8"))
ex.clusters <- ex.se.obj@ident; names(ex.clusters) <- ex.se.obj@cell.names; levels(ex.clusters);
ex.clusters <- plyr::revalue(ex.clusters, replace = 
                               c("Ex1" = "L2/3", "Ex3a" = "L4", "Ex3b" = "L4", "Ex3c" = "L4", "Ex3d" = "L4",
                                 "Ex4" = "L4/5", "Ex5a" = "L5", "Ex5b" = "L5", "Ex6a" = "L6", "Ex6b" = "L6", 
                                 "Ex8" = "L6b"))
levels(ex.clusters)

ex.se.obj <- FindVariableGenes(ex.se.obj, x.low.cutoff = 0.025, y.cutoff = 0.5)
length(ex.se.obj@var.genes)

ex.se.obj <- ScaleData(ex.se.obj, genes.use = ex.se.obj@var.genes, vars.to.regress = "nUMI", model.use = "linear")
ex.se.obj <- RunPCA(ex.se.obj, pcs.compute = 40, do.print = F)
PCElbowPlot(ex.se.obj, num.pc = 40)

pcs.use <- 20
ex.se.obj <- RunTSNE(ex.se.obj, dims.use = 1:pcs.use)

## Run NMF
ex.norm.counts <- ScaleCounts(se.obj@raw.data[,ex.se.obj@cell.names], method = "ft")
ex.var.genes <- intersect(ex.se.obj@var.genes, rownames(ex.norm.counts))

ex.k.res <- FindNumFactors(ex.norm.counts[ex.var.genes,], k.range = k.range, n.cores = n.cores, 
                           na.frac = 0.25, seed = seed, loss = loss)
ex.k.res$k

k <- 14
ex.nmf.res <- RunNMF(ex.norm.counts[ex.var.genes,], k = k, alpha = 0, init = "ica", n.cores = n.cores,
                     loss = loss, init.zeros = "random")
ex.nmf.res$W <- ProjectFeatures(ex.norm.counts, ex.nmf.res$H, loss = loss, n.cores = n.cores)
ex.nmf.scores <- ex.nmf.res$H

ex.se.obj <- BuildSNN(ex.se.obj, k = 20, prune.SNN = 1/20, dims.use = 1:pcs.use, force.recalc = T)
ex.swne.embedding <- EmbedSWNE(ex.nmf.scores, ex.se.obj@snn, alpha.exp = 5.0, snn.exp = 0.25, 
                               n_pull = 4, snn.factor.proj = T, pca.red = F, dist.use = "cosine")

## Hide all factors
ex.swne.embedding$H.coords$name <- ""

genes.embed <- c("NTNG1", "DAB1", "HS3ST2", "DCC", "POSTN")
ex.swne.embedding <- EmbedFeatures(ex.swne.embedding, ex.nmf.res$W, genes.embed, n_pull = 4, scale.cols = F)

pdf("snDropSeq_swne_ex_neurons.pdf", width = 6.5, height = 6.5)
PlotSWNE(ex.swne.embedding, alpha.plot = 0.3, sample.groups = ex.clusters, do.label = T, 
         label.size = 4.0, pt.size = 1.0, show.legend = F, seed = seed)
dev.off()

pdf("snDropSeq_swne_ex_neurons_nolabel.pdf", width = 6.5, height = 6.5)
PlotSWNE(ex.swne.embedding, alpha.plot = 0.3, sample.groups = ex.clusters, do.label = T, 
         label.size = 0, pt.size = 1.0, show.legend = F, seed = seed)
dev.off()

ex.tsne.scores <- GetCellEmbeddings(ex.se.obj, reduction.type = "tsne")
pdf("snDropSeq_tsne_ex_neurons.pdf", width = 6, height = 6)
PlotDims(ex.tsne.scores, sample.groups = ex.clusters, pt.size = 0.5, label.size = 4, 
         alpha.plot = 0.35, show.legend = F, show.axes = F, seed = seed)
dev.off()

pdf("snDropSeq_tsne_ex_neurons_nolabel.pdf", width = 6, height = 6)
PlotDims(ex.tsne.scores, sample.groups = ex.clusters, pt.size = 0.5, label.size = 0, 
         alpha.plot = 0.4, show.legend = F, show.axes = F, seed = seed)
dev.off()


## Find layer specific markers
ex.se.obj <- SetIdent(ex.se.obj, ident.use = ex.clusters)
ex.layer.markers.df <- FindAllMarkers(ex.se.obj, logfc.threshold = 0.02, only.pos = T)
ex.top.markers.df <- Reduce("rbind", by(ex.layer.markers.df, ex.layer.markers.df$cluster, head, n = 3))

ex.layer.markers.mat <- UnflattenDataframe(ex.layer.markers.df, output.name = "avg_logFC",
                                           row.col = "gene", col.col = "cluster")
ex.layer.markers.mat[is.na(ex.layer.markers.mat)] <- 0
ex.layer.markers.mat <- ex.layer.markers.mat[unique(ex.top.markers.df$gene),]

pdf("snDropSeq_ex_gene_layer_heatmap.pdf", width = 4.5, height = 4.5)
ggHeat(ex.layer.markers.mat, clustering = "both", x.lab.size = 14, y.lab.size = 12)
dev.off()


ex.top.genes.df <- SummarizeAssocFeatures(ex.gene.factor.assoc, features.return = 6)
ex.top.genes.df <- subset(ex.top.genes.df, assoc_score > 0.3 & 
                            factor %in% c("factor_6", "factor_12"))
ex.gene.nmf.heat <- ex.gene.factor.assoc[unique(ex.top.genes.df$feature),]
pdf("snDropSeq_ex_metagene_heatmap.pdf", width = 5, height = 3.5)
ggHeat(ex.gene.nmf.heat, clustering = "both", x.lab.size = 11, y.lab.size = 10)
dev.off()


## Validate gene embeddings
gene <- "DAB1"
ex.gene.swne.embedding <- ex.swne.embedding
ex.gene.swne.embedding$H.coords$name <- ""
ex.gene.swne.embedding$feature.coords <- subset(ex.gene.swne.embedding$feature.coords, name == gene)

pdf("snDropSeq_ex_DAB1_feature_plot.pdf", width = 4.5, height = 4.5)
FeaturePlotSWNE(ex.gene.swne.embedding, norm.counts[gene,], pt.size = 0.5, label.size = 0)
dev.off()

save.image("snDropSeq_swne.RData")
