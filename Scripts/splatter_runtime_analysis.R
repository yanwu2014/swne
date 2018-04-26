## Test SWNE runtime
library(splatter)
library(pagoda2)
setwd("/media/Scratch_SSD/Yan/R_Analysis/swne-dev/simulations")

## Generate data
splat.seed <- 3250879
n.cells <- c(1e3, 5e3, 1e4, 5e4)

simulated.data <- lapply(n.cells, function(n) {
  splat.params <- newSplatParams(batchCells = n, nGenes = 1e4,
                                 group.prob = c(0.45, 0.15, 0.15, 0.15, 0.1),
                                 de.prob = c(0.25, 0.05, 0.05, 0.05, 0.25),
                                 seed = splat.seed)

  ## Generate synthetic data using fitted parameters
  splat.groups <- splatSimulateGroups(splat.params, verbose = F)
  counts <- as(splat.groups@assays@.xData$data$counts, "dgCMatrix")
  colnames(counts) <- paste(colnames(counts), splat.groups@colData@listData$Group, sep = "_")

  return(as(counts, "dgCMatrix"))
})

library(swne)
run_swne <- function(counts, n.genes, init, k, n.cores) {
  counts <- ScaleCounts(counts, method = "ft")
  
  varinfo <- AdjustVariance(counts, gam.k = 10)
  varinfo <- varinfo[order(varinfo$lp),]
  od.genes <- rownames(varinfo[1:n.genes,])
  
  nmf <- RunNMF(counts[od.genes,], init = init, k = k, n.cores = n.cores)
  
  pca.center <- Matrix::rowMeans(counts[od.genes,])
  pca <- irlba::irlba(Matrix::t(counts[od.genes,]) - pca.center, nv = 40, maxit = 500)
  pca.red <- t(pca$u); colnames(pca.red) <- colnames(counts);
  
  snn <- CalcSNN(pca.red, k = 30)
  
  emb <- EmbedSWNE(nmf$H, snn, n_pull = 4)
  return(emb)
}

library(compiler)
run_swne_cmp <- cmpfun(run_swne)

n.features <- 3e3
init <- "ica"
k <- 16
n.cores <- 16

library(microbenchmark)
runtimes.ncells <- lapply(simulated.data, function(counts) {
  microbenchmark(run_swne_cmp(counts, n.features, init, k = k, n.cores = n.cores), times = 4)
})

ngenes.range <- c(3e3, 1e3, 500, 200)
runtimes.ngenes <- lapply(ngenes.range, function(n) {
  microbenchmark(run_swne_cmp(simulated.data[[3]], n, init, k = k, n.cores = n.cores), times = 4)
})

runtimes.ncells.df <- data.frame(mean = sapply(runtimes.ncells, function(df) mean(df$time)/1e9),
                                 sd = sapply(runtimes.ncells, function(df) sd(df$time)/1e9))
runtimes.ncells.df$cells <- n.cells

runtimes.ngenes.df <- data.frame(mean = sapply(runtimes.ngenes, function(df) mean(df$time)/1e9),
                                 sd = sapply(runtimes.ngenes, function(df) sd(df$time)/1e9))
runtimes.ngenes.df$genes <- ngenes.range

library(ggplot2)

pdf("splatter_runtime_ncells.pdf", width = 4, height = 4)
ggplot(runtimes.ncells.df, aes(x = cells, y = mean)) + 
  geom_line() + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.05,
                position = position_dodge()) +
  scale_y_log10(breaks = c(50, 200, 500, 1000)) + scale_x_log10(breaks = n.cells) + 
  theme_classic() + xlab("Number of cells") + ylab("Runtime (sec)") + 
  theme(axis.text = element_text(size = 11))
dev.off()

pdf("splatter_runtime_ngenes.pdf", width = 4, height = 4)
ggplot(runtimes.ngenes.df, aes(x = genes, y = mean)) + 
  geom_line() + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.05,
                position = position_dodge()) +
  scale_y_log10(breaks = c(25, 50, 100, 200)) + scale_x_log10(breaks = ngenes.range) + 
  theme_classic() + xlab("Number of genes") + ylab("Runtime (sec)") +
  theme(axis.text = element_text(size = 11))
dev.off()

save.image("splatter.runtime.sim.RData")
