library(Seurat)
library(swne)

# sim.groups <- splatter::splatSimulateGroups(nGenes = 6000, groupCells = c(200, 200, 200, 400),
#                                             de.prob = c(0.025, 0.025, 0.025, 0.5), verbose = F)
# sim.counts <- as.matrix(sim.groups@assayData$counts)
# sim.metadata <- sim.groups@phenoData@data
# 
# write.table(as.matrix(sim.counts), file = "sim.counts.tsv", sep = "\t")
# write.table(sim.metadata, file = "sim.metadata.tsv", sep = "\t")

n.cores <- 16

sim.counts <- as(as.matrix(read.table("sim.counts.tsv", sep = "\t", header = T)), "dgCMatrix")
sim.metadata <- read.table("sim.metadata.tsv", sep = "\t", header = T)

se.obj <- CreateSeuratObject(sim.counts, meta.data = sim.metadata)
se.obj <- SetAllIdent(se.obj, id = "Group")
se.obj <- NormalizeData(se.obj, scale.factor = median(Matrix::colSums(sim.counts)))
se.obj <- ScaleData(se.obj, model.use = "linear", genes.use = se.obj@var.genes)

se.obj <- RunPCA(se.obj, pc.genes = rownames(se.obj@scale.data), do.print = F)
PCElbowPlot(se.obj, num.pc = 20)

dims.use <- 5
se.obj <- BuildSNN(se.obj, dims.use = 1:dims.use, k.param = 10, k.scale = 10, force.recalc = T)
se.obj <- RunTSNE(se.obj, dims.use = 1:dims.use)

sim.norm.counts <- ScaleCounts(sim.counts, batch = NULL, method = "ft", adj.var = T)
cell.clusters <- factor(sim.metadata$Group); names(cell.clusters) <- colnames(sim.counts);

nmf.res <- RunNMF(sim.norm.counts[se.obj@var.genes,], k = dims.use, alpha = 0, init = "ica", 
                  n.cores = n.cores, loss = "mse", init.zeros = "random")
nmf.res$W <- ProjectFeatures(sim.norm.counts, nmf.res$H, loss = "mse", n.cores = n.cores)
H <- nmf.res$H

top.assoc.genes.df <- SummarizeAssocFeatures(nmf.res$W, features.return = 2)

swne.embedding <- EmbedSWNE(H, SNN = se.obj@snn, alpha.exp = 1, snn.exp = 1, n_pull = 4, dist.use = "IC")
swne.embedding <- RenameFactors(swne.embedding, c("factor_1" = "biological_factor_1", "factor_2" = "biological_factor_2"))
swne.embedding <- EmbedFeatures(swne.embedding, nmf.res$W, c("Gene5707", "Gene4298"), n_pull = 4)

seed <- 43859279
pdf("sim_swne_plot_nolabels.pdf", width = 5, height = 5)
PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = cell.clusters, do.label = T, label.size = 0, 
         pt.size = 1.5, seed = seed, show.legend = F)
dev.off()

pc.scores <- GetCellEmbeddings(se.obj, reduction.type = "pca", dims.use = 1:2)
pdf("sim_pca_plot_nolabels.pdf", width = 5, height = 5)
PlotDims(pc.scores, sample.groups = cell.clusters, x.lab = NULL, y.lab = NULL, seed = seed, 
         show.axes = T, alpha.plot = 0.5, label.size = 0, show.legend = F) + 
  theme(axis.ticks = element_blank(), axis.text = element_blank())
dev.off()

tsne.scores <- GetCellEmbeddings(se.obj, reduction.type = "tsne", dims.use = 1:2)
pdf("sim_tsne_plot_nolabels.pdf", width = 5, height = 5)
PlotDims(tsne.scores, sample.groups = cell.clusters, x.lab = NULL, y.lab = NULL, seed = seed, 
         show.axes = T, alpha.plot = 0.5, label.size = 0, show.legend = F) +
  theme(axis.ticks = element_blank(), axis.text = element_blank())
dev.off()

save.image("sim_swne.RData")
