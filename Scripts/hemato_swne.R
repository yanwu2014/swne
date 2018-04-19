## Read data
library(Seurat)
library(swne)

## Create Seurat object and create SNN
counts <- ReadData("hemato_expr_debatched.tsv", make.sparse = T)
info.genes <- scan("hemato_info_genes.txt", sep = "\n", what = character())

clusters.df <- read.table("hemato_cluster_mapping.csv", sep = ",")
clusters <- clusters.df[[2]]; names(clusters) <- clusters.df[[1]];

counts <- counts[,names(clusters)]
counts <- FilterData(counts, min.samples.frac = 0.005, trim = 0.005, min.nonzero.features = 300)
info.genes <- info.genes[info.genes %in% rownames(counts)]
dim(counts)

se.obj <- CreateSeuratObject(counts)
se.obj <- SetIdent(se.obj, cells.use = names(clusters), clusters)
se.obj@ident <- plyr::revalue(se.obj@ident, 
                              c("1" = 'Ery', "2" = 'Ery', "3" = 'Ery', "4" = 'Ery', "5" = 'Ery',
                                "6" = 'Ery', "7" = 'MP/EP', "8" = 'MK', "9" = 'GMP', "10" = 'GMP',
                                "11" = 'DC', "12" = 'Bas', "13" = 'Bas', "14" = 'M', "15" = 'M',
                                "16" = 'Neu', "17" = 'Neu', "18" = 'Eos', "19" = 'lymphoid'))
levels(se.obj@ident)
rm(counts)

se.obj <- SubsetData(se.obj, ident.remove = "lymphoid")
se.obj <- NormalizeData(se.obj, normalization.method = "LogNormalize")
se.obj@var.genes <- info.genes

se.obj <- ScaleData(se.obj, genes.use = se.obj@var.genes, model.use = "negbinom")
se.obj <- RunPCA(se.obj, pc.genes = rownames(counts), do.print = F, pcs.compute = 40)
PCElbowPlot(se.obj, num.pc = 40)

se.obj <- RunTSNE(se.obj, dims.use = 1:16)
se.obj <- BuildSNN(se.obj, dims.use = 1:16, k.param = 30, prune.SNN = 1/30, force.recalc = T)
clusters <- se.obj@ident; names(clusters) <- se.obj@cell.names;

## SWNE parameters
k.range <- seq(2,20,2)
n.cores <- 16

norm.counts <- ScaleCounts(se.obj@raw.data, method = "ft")[,names(clusters)]
k.res <- FindNumFactors(norm.counts[se.obj@var.genes,], k.range = k.range, n.cores = n.cores, na.frac = 0.2,
                        seed = 32590, do.plot = T, loss = loss, recon.err = loss)
k.res$k

k <- 14
nmf.res <- RunNMF(norm.counts[se.obj@var.genes,], k = k, alpha = 0, init = "ica", n.cores = n.cores,
                  loss = loss, init.zeros = "random")
nmf.res$W <- ProjectFeatures(norm.counts, nmf.res$H, loss = "mse", n.cores = n.cores)
nmf.scores <- nmf.res$H

snn.matrix <- se.obj@snn
swne.embedding <- EmbedSWNE(nmf.scores, SNN = snn.matrix, alpha.exp = 2.0, snn.exp = 0.1, n_pull = 3, 
                            dist.use = "cosine")

swne.embedding <- RenameFactors(swne.embedding, name.mapping = c("factor_4" = "Metal binding",
                                                                 "factor_8" = "Antigen presentation",
                                                                 "factor_10" = "Platelet generation"))

## Embed selected genes onto swne plot
genes.embed <- c("Apoe", "Mt2","Flt3", "Sun2", "Pglyrp1")
swne.embedding <- EmbedFeatures(swne.embedding, nmf.res$W, genes.embed, n_pull = 4, scale.cols = F)

cluster_colors <- c("Bas" = "#ff6347", "DC" = "#46C7EF", "Eos" = "#EFAD1E", "Ery" = "#8CB3DF", "M" = "#53C0AD",
                    "MP/EP" = "#4EB859", "GMP" = "#D097C4", "MK" = "#ACC436", "Neu" = "#F5918A")

pdf("hemato_swne_plot.pdf", width = 6, height = 6)
PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = clusters, do.label = T,
         label.size = 3.5, pt.size = 1.5, show.legend = F) + 
  scale_color_manual(values = cluster_colors)
dev.off()

pdf("hemato_swne_plot_nolabels.pdf", width = 6, height = 6)
PlotSWNE(swne.embedding, alpha.plot = 0.5, sample.groups = clusters, do.label = F,
         label.size = 0, pt.size = 1.5, show.legend = F) + 
  scale_color_manual(values = cluster_colors)
dev.off()

tsne.scores <- GetCellEmbeddings(se.obj, reduction.type = "tsne")
pdf("hemato_tsne_plot.pdf", width = 6, height = 6)
PlotDims(tsne.scores, sample.groups = clusters, show.legend = F, show.axes = F, 
         alpha.plot = 0.5, label.size = 3.5, pt.size = 1.5) + 
  scale_color_manual(values = cluster_colors)
dev.off()

pdf("hemato_swne_plot_nolabels.pdf", width = 6, height = 6)
PlotSWNE(swne.embedding, alpha.plot = 0.5, sample.groups = clusters, do.label = F,
         label.size = 0, pt.size = 1.5, show.legend = F) + 
  scale_color_manual(values = cluster_colors)
dev.off()

pdf("hemato_tsne_plot_nolabels.pdf", width = 6, height = 6)
PlotDims(tsne.scores, sample.groups = clusters, show.legend = F, show.axes = F, 
         alpha.plot = 0.5, pt.size = 1.5, label.size = 0) + 
  scale_color_manual(values = cluster_colors)
dev.off()

## Overlay SWNE plot with pseudotime
ps.time <- read.table("hemato_monocle_pseudotime.tsv", sep = "\t", header = T, stringsAsFactors = F)
pseudotime <- ps.time$pseudotime; names(pseudotime) <- ps.time$cell.names;

pdf("hemato_swne_pseudotime_plot.pdf", width = 6.5, height = 6)
FeaturePlotSWNE(swne.embedding, pseudotime, alpha.plot = 0.4, label.size = 3.5, pt.size = 1.5)
dev.off()

pdf("hemato_tsne_pseudotime_plot.pdf", width = 6.5, height = 6)
FeaturePlotDims(tsne.scores, pseudotime, alpha.plot = 0.4, show.axes = F, pt.size = 1.5)
dev.off()

pdf("hemato_swne_pseudotime_plot_nolabels.pdf", width = 6.5, height = 6)
FeaturePlotSWNE(swne.embedding, pseudotime, alpha.plot = 0.4, label.size = 0, pt.size = 1.5)
dev.off()


## Associate factors with cell clusters
clusters.list <- UnflattenGroups(clusters)
clusters.matrix <- t(swne:::.genesets_indicator(clusters.list, inv = F, return.numeric = T))
cluster.nmf.assoc <- FactorAssociation(clusters.matrix, nmf.scores, n.cores = n.cores, metric = "IC")

pdf("hemato_cluster_factor_heatmap.pdf", width = 6.0, height = 4.0)
ggHeat(cluster.nmf.assoc, clustering = "both", x.lab.size = 14, y.lab.size = 12)
dev.off()

## Associate factors with genes using the gene loadings (W) matrix
# gene.loadings <- t(apply(nmf.res$W, 1, function(x) (x - min(x))/(max(x) - min(x))))
gene.loadings <- nmf.res$W
top.factor.genes.df <- SummarizeAssocFeatures(gene.loadings, features.return = 8)
write.table(top.factor.genes.df, file = "hemato_factor_markers.tsv", sep = "\t")

## Make gene factor association heatmaps
gene.factor.heat <- gene.loadings[unique(top.factor.genes.df$feature),]

pdf("hemato_gene_factor_heatmap.pdf", width = 6.0, height = 6.25)
ggHeat(gene.factor.heat, clustering = "both", x.lab.size = 14, y.lab.size = 12)
dev.off()

cluster.genes.df <- FindAllMarkers(se.obj, logfc.threshold = 0.0, only.pos = T)
cluster.genes.mat <- UnflattenDataframe(cluster.genes.df, output.name = "avg_logFC",
                                        row.col = "gene", col.col = "cluster")
cluster.genes.mat[is.na(cluster.genes.mat)] <- 0

top.cluster.genes.df <- Reduce("rbind", by(cluster.genes.df, cluster.genes.df$cluster, head, n = 3))
cluster.genes.heat <- cluster.genes.mat[unique(top.cluster.genes.df$gene),]

pdf("hemato_gene_cluster_heatmap.pdf", width = 6.0, height = 6)
ggHeat(cluster.genes.heat, clustering = "both", x.lab.size = 14, y.lab.size = 12)
dev.off()


## Validate gene embeddings
gene <- "Sun2"
gene.swne.embedding <- swne.embedding
gene.swne.embedding$H.coords$name <- ""
gene.swne.embedding$feature.coords <- subset(gene.swne.embedding$feature.coords, name == gene)

pdf("hemato_Sun2_feature_plot.pdf", width = 4.5, height = 4.5)
FeaturePlotSWNE(gene.swne.embedding, norm.counts[gene,], label.size = 0)
dev.off()

save.image("hemato_swne.RData")
