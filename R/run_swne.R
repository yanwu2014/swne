#' Wrapper for the running SWNE analysis functions
#'
#' @param object A seurat-class object (normalised)
#' @param dist.use Similarity function to use for calculating factor positions (passed to EmbedSWNE). Options include pearson (correlation), IC (mutual information), cosine, euclidean.
#' @param n.cores Number of cores to use (passed to FindNumFactors)
#' @param k Number of NMF factors (passed to RunNMF). If none given, will be derived from FindNumFactors.
#'
#' @return A list of factor (H.coords) and sample coordinates (sample.coords) in 2D
#'
#' @export
RunSWNE <- function(object, dist.metric = "euclidean", n.cores = 3, k){
  #Run SWNE with correlation matrix
  object_norm <- ExtractNormCounts(object, obj.type = "seurat", rescale = F, rescale.method = "log", batch = NULL)
  var_genes <- intersect(object@var.genes, rownames(object_norm));
  length(var_genes)
  cell_clusters <- object@ident
  names(cell_clusters) <- object@cell.names
  levels(cell_clusters)
  loss <- "mse" ## Loss function
  if(missing(k)){
    n.cores <- n.cores ## Number of cores to use
    k.range <- seq(2,10,2) ## Range of factors to iterate over
    k.res <- FindNumFactors(object_norm[var_genes,], k.range = k.range, n.cores = n.cores, do.plot = F, loss = loss)
    print(k.res$k, "factors")
    k <- k.res$k
  }
  if(k < 3) warning("k must be an integer of 3 or higher")
  k <- max(k, 3)
  nmf.res <- RunNMF(object_norm[var_genes,], k = k, alpha = 0, init = "ica", n.cores = n.cores, loss = loss)
  nmf.scores <- nmf.res$H
  # pc.scores <- t(GetCellEmbeddings(test, reduction.type = "pca", dims.use = 1:k))
  # snn <- CalcSNN(pc.scores)
  snn <- object@snn
  #correct for aggregrated cell barcodes
  colnames(nmf.scores) <- colnames(snn)
  alpha.exp <- 1.25 # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
  snn.exp <- 1.0 # Lower this < 1.0 to move similar cells closer to each other
  n_pull <- max(k.res$k, 4) # The number of factors pulling on each cell. Must be at least 3.
  swne_embedding <- EmbedSWNE(nmf.scores, snn, alpha.exp = alpha.exp, snn.exp = snn.exp,
                              n_pull = n_pull, dist.use = dist.metric)
  return(swne_embedding)
}
