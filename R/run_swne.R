#' Wrapper for the running SWNE analysis functions
#'
#' @param object A seurat-class object (normalised)
#' @param data.matrix a data matrix (genes x cells) which has been pre-normalised
#' @param dist.use Similarity function to use for calculating factor positions (passed to EmbedSWNE). Options include pearson (correlation), IC (mutual information), cosine, euclidean.
#' @param n.cores Number of cores to use (passed to FindNumFactors)
#' @param k Number of NMF factors (passed to RunNMF). If none given, will be derived from FindNumFactors.
#' @param var.genes vector to specify variable genes. Will infer from Suerat or use full dataset if not given.
#' @param loss loss function to use (passed to RunNMF)
#' @param alpha.exp Increasing alpha.exp increases how much the NMF factors "pull" the samples (passed to EmbedSWNE)
#' @param snn.exp Decreasing snn.exp increases the effect of the similarity matrix on the embedding (passed to EmbedSWNE)
#' @param var.exp Proportion of genes selected from most variable
#'
#' @return A list of factor (H.coords) and sample coordinates (sample.coords) in 2D
#'
#' @export RunSWNE
#' @rdname RunSWNE
#' @usage NULL
RunSWNE <- function(x, ...) {
  UseMethod("RunSWNE")
}


#' @rdname RunSWNE
#' @method RunSWNE seurat
#' @export
#'
RunSWNE.seurat <- function(object, dist.metric = "euclidean", n.cores = 3, k, var.genes, loss = "mse",
                           alpha.exp = 1.25, # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
                           snn.exp = 1.0 # Lower this < 1.0 to move similar cells closer to each other
){
  object_norm <- ExtractNormCounts(object, obj.type = "seurat", rescale = F, rescale.method = "log", batch = NULL)
  if(missing(var.genes)) var_genes <- intersect(object@var.genes, rownames(object_norm))
  print(paste(length(var_genes), "variable genes"))
  cell_clusters <- object@ident
  names(cell_clusters) <- object@cell.names
  print(paste(levels(cell_clusters), "cell clusters"))
  if(missing(k)){
    n.cores <- n.cores ## Number of cores to use
    k.range <- seq(2,10,2) ## Range of factors to iterate over
    k.res <- FindNumFactors(object_norm[var_genes,], k.range = k.range, n.cores = n.cores, do.plot = F, loss = loss)
    print(paste(k.res$k, "factors"))
    k <- k.res$k
  }
  if(k < 3) warning("k must be an integer of 3 or higher")
  k <- max(k, 3)
  nmf.res <- RunNMF(object_norm[var_genes,], k = k, alpha = 0, init = "ica", n.cores = n.cores, loss = loss)
  nmf.scores <- nmf.res$H
  if(sum(dim(object@snn)) < 2){
    object <- RunPCA(object)
    object <- FindClusters(object = object, reduction.type = "pca", dims.use = 1:10,
                           resolution = 0.6, print.output = 0, save.SNN = TRUE)
    if(sum(dim(object@snn)) < 2){
      pc.scores <- t(GetCellEmbeddings(object, reduction.type = "pca", dims.use = 1:k))
      object@snn <- CalcSNN(pc.scores)
    }
  }
  snn <- object@snn
  #correct for aggregrated cell barcodes
  colnames(nmf.scores) <- colnames(snn)
  alpha.exp <- 1.25 # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
  snn.exp <- 1.0 # Lower this < 1.0 to move similar cells closer to each other
  n_pull <- k # The number of factors pulling on each cell. Must be at least 3.
  swne_embedding <- EmbedSWNE(nmf.scores, snn, alpha.exp = alpha.exp, snn.exp = snn.exp,
                              n_pull = n_pull, dist.use = dist.metric)
  return(swne_embedding)
}

#' @rdname RunSWNE
#' @method RunSWNE matrix
#' @export
RunSWNE.matrix <- function(data.matrix, dist.metric = "euclidean", n.cores = 3, k, var.genes = rownames(data.matrix), loss = "mse",
                            alpha.exp = 1.25, # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
                            snn.exp = 1.0, # Lower this < 1.0 to move similar cells closer to each other
                            var.exp = 0.05 # Increase this to include more genes
){
  object_norm <- data.matrix - min(data.matrix)
  vars <- apply(data.matrix, 1, var)
  var_genes <- rownames(object_norm)[vars >= quantile(vars, 1 - var.exp)]
  print(paste(length(var_genes), "variable genes"))
  if(missing(k)){
    n.cores <- n.cores ## Number of cores to use
    k.range <- seq(2,10,2) ## Range of factors to iterate over
    k.res <- FindNumFactors(object_norm[var_genes,], k.range = k.range, n.cores = n.cores, do.plot = F, loss = loss)
    print(paste(k.res$k, "factors"))
    k <- k.res$k
  }
  if(k < 3) warning("k must be an integer of 3 or higher")
  k <- max(k, 3)
  nmf.res <- RunNMF(object_norm[var_genes,], k = k, alpha = 0, init = "ica", n.cores = n.cores, loss = loss)
  nmf.scores <- nmf.res$H
  pc.scores <- prcomp(data.matrix, center = TRUE, scale = TRUE)$rotation[,1:k]
  snn <- CalcSNN(t(pc.scores), k = k)
  #correct for aggregrated cell barcodes
  colnames(nmf.scores) <- rownames(snn) <- colnames(snn)
  n_pull <- k # The number of factors pulling on each cell. Must be at least 3.
  swne_embedding <- EmbedSWNE(nmf.scores, snn, alpha.exp = alpha.exp, snn.exp = snn.exp,
                              n_pull = k, dist.use = dist.metric)
  return(swne_embedding)
}

#' @rdname RunSWNE
#' @method RunSWNE default
#' @export
RunSWNE.default <- function(data.matrix, dist.metric = "euclidean", n.cores = 3, k, var.genes = rownames(data.matrix), loss = "mse",
                            alpha.exp = 1.25, # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
                            snn.exp = 1.0 # Lower this < 1.0 to move similar cells closer to each other
){
  object_norm <- data.matrix
  var_genes <- intersect(var.genes, rownames(object_norm))
  print(paste(length(var_genes), "variable genes"))
  if(missing(k)){
    n.cores <- n.cores ## Number of cores to use
    k.range <- seq(2,10,2) ## Range of factors to iterate over
    k.res <- FindNumFactors(object_norm[var_genes,], k.range = k.range, n.cores = n.cores, do.plot = F, loss = loss)
    print(paste(k.res$k, "factors"))
    k <- k.res$k
  }
  if(k < 3) warning("k must be an integer of 3 or higher")
  k <- max(k, 3)
  nmf.res <- RunNMF(object_norm[var_genes,], k = k, alpha = 0, init = "ica", n.cores = n.cores, loss = loss)
  nmf.scores <- nmf.res$H
  pc.scores <- prcomp(data.matrix, center = TRUE, scale = TRUE)$rotation[,1:k]
  snn <- CalcSNN(t(pc.scores), k = k)
  #correct for aggregrated cell barcodes
  colnames(nmf.scores) <- rownames(snn) <- colnames(snn)
  n_pull <- k # The number of factors pulling on each cell. Must be at least 3.
  swne_embedding <- EmbedSWNE(nmf.scores, snn, alpha.exp = alpha.exp, snn.exp = snn.exp,
                              n_pull = k, dist.use = dist.metric)
  return(swne_embedding)
}
