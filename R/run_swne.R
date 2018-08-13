#Wrapper for the "swne" package
## object (seurat object)
## dist.metric distance to pass to EmbedSWNE
## n.cores number of cores to use for FindNumFactors
#object seurat-class object
RunSWNE <- function(object, dist.metric = "euclidean", n.cores = 3){
  #Run SWNE with correlation matrix
  object_norm <- ExtractNormCounts(object, obj.type = "seurat", rescale = F, rescale.method = "log", batch = NULL)
  var_genes <- intersect(object@var.genes, rownames(object_norm));
  length(var_genes)
  cell_clusters <- object@ident
  names(cell_clusters) <- object@cell.names
  levels(cell_clusters)
  loss <- "mse" ## Loss function
  n.cores <- n.cores ## Number of cores to use
  k.range <- seq(2,10,2) ## Range of factors to iterate over
  k.res <- FindNumFactors(object_norm[var_genes,], k.range = k.range, n.cores = n.cores, do.plot = T, loss = loss)
  print(k.res$k, "factors")
  k <- max(k.res$k, 3)
  nmf.res <- RunNMF(object_norm[var_genes,], k = k, alpha = 0, init = "ica", n.cores = n.cores, loss = loss)
  nmf.scores <- nmf.res$H
  # pc.scores <- t(GetCellEmbeddings(test, reduction.type = "pca", dims.use = 1:k))
  # snn <- CalcSNN(pc.scores)
  snn <- object@snn
  alpha.exp <- 1.25 # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
  snn.exp <- 1.0 # Lower this < 1.0 to move similar cells closer to each other
  n_pull <- max(k.res$k, 4) # The number of factors pulling on each cell. Must be at least 3.
  swne_embedding <- EmbedSWNE(nmf.scores, snn, alpha.exp = alpha.exp, snn.exp = snn.exp,
                              n_pull = n_pull, dist.use = dist.metric)
}
