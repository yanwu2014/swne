#' Wrapper for running SWNE analysis
#'
#' @param object A Seurat or Pagoda2 object with normalized data
#' @param data.matrix a data matrix (genes x cells) which has been pre-normalized
#' @param batch Vector of batch effects to correct for
#' @param proj.method Method to use to project factors in 2D. Either "sammon" or "umap"
#' @param dist.use Similarity function to use for calculating factor positions (passed to EmbedSWNE).
#'                 Options include pearson (correlation), IC (mutual information), cosine, euclidean.
#' @param n.cores Number of cores to use (passed to FindNumFactors)
#' @param k Number of NMF factors (passed to RunNMF). If none given, will be derived from FindNumFactors.
#' @param k.range Range of factors for FindNumFactors to iterate over if k is not given
#' @param var.genes vector to specify variable genes. Will infer from Seurat or use full dataset if not given.
#' @param loss loss function to use (passed to RunNMF)
#' @param alpha.exp Increasing alpha.exp increases how much the NMF factors "pull" the samples (passed to EmbedSWNE)
#' @param snn.exp Decreasing snn.exp increases the effect of the similarity matrix on the embedding (passed to EmbedSWNE)
#' @param n.var.genes Number of variable genes to use
#' @param n_pull Maximum number of factors "pulling" on each sample
#' @param genes.embed Genes to add to the SWNE embedding
#' @param hide.factors Hide factors when plotting SWNE embedding
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
RunSWNE.seurat <- function(object, proj.method = "umap", dist.metric = "cosine", n.cores = 8, k, k.range, var.genes,
                           loss = "mse", genes.embed, hide.factors = T, n_pull = 3,
                           alpha.exp = 1.25, # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
                           snn.exp = 1.0 # Lower this < 1.0 to move similar cells closer to each other
){
  object_norm <- ExtractNormCounts(object, obj.type = "seurat", rescale = F, rescale.method = "log", batch = NULL)

  if (missing(var.genes)) var.genes <- intersect(object@var.genes, rownames(object_norm))
  var.genes <- intersect(var.genes, rownames(object_norm))
  print(paste(length(var.genes), "variable genes to use"))

  if (missing(k)) {
    if (missing(k.range)) k.range <- seq(2,20,2) ## Range of factors to iterate over
    k.res <- FindNumFactors(object_norm[var.genes,], k.range = k.range, n.cores = n.cores, do.plot = F, loss = loss)
    print(paste(k.res$k, "factors")); k <- k.res$k;
  }

  if (k < 3) {
    warning("k must be an integer of 3 or higher")
    k <- 3
  }

  if(sum(dim(object@snn)) < 2){
    object <- RunPCA(object, pc.genes = var.genes, do.print = F, pcs.compute = max(k,20))
    pc.scores <- t(GetCellEmbeddings(object, reduction.type = "pca", dims.use = 1:k))
    snn <- CalcSNN(pc.scores, k = 20, prune.SNN = 1/20)
  } else {
    snn <- object@snn
  }

  if (missing(genes.embed)) genes.embed <- NULL
  run_swne(object_norm, var.genes, snn, k, alpha.exp, snn.exp, n_pull, proj.method, dist.metric, genes.embed,
           loss, n.cores, hide.factors)
}



#' @rdname RunSWNE
#' @method RunSWNE Pagoda2
#' @export
#'
RunSWNE.Pagoda2 <- function(object, proj.method = "umap", dist.metric = "cosine", n.cores = 8, k, k.range, var.genes,
                            loss = "mse", genes.embed, hide.factors = T, n_pull = 3, n.var.genes = 3000,
                            alpha.exp = 1.25, # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
                            snn.exp = 1.0 # Lower this < 1.0 to move similar cells closer to each other
){
  object_norm <- ExtractNormCounts(object, obj.type = "pagoda2", rescale = F, rescale.method = "log", batch = NULL)

  if (missing(var.genes)) var.genes <- rownames(p2$misc$varinfo[order(p2$misc$varinfo$lp),])[1:n.var.genes]
  var.genes <- intersect(var.genes, rownames(object_norm))
  print(paste(length(var.genes), "variable genes to use"))

  if (missing(k)) {
    if (missing(k.range)) k.range <- seq(2,20,2) ## Range of factors to iterate over
    k.res <- FindNumFactors(object_norm[var.genes,], k.range = k.range, n.cores = n.cores, do.plot = F, loss = loss)
    print(paste(k.res$k, "factors")); k <- k.res$k;
  }

  if (k < 3) {
    warning("k must be an integer of 3 or higher")
    k <- 3
  }

  object$calculatePcaReduction(nPcs = max(k,20), odgenes = var.genes)
  pc.scores <- t(object$reductions$PCA[,1:k])
  snn <- CalcSNN(pc.scores, k = 20, prune.SNN = 1/20)

  if (missing(genes.embed)) genes.embed <- NULL
  run_swne(object_norm, var.genes, snn, k, alpha.exp, snn.exp, n_pull, proj.method, dist.metric, genes.embed,
           loss, n.cores, hide.factors)
}



#' @rdname RunSWNE
#' @method RunSWNE dgCMatrix
#' @export
RunSWNE.dgCMatrix <- function(data.matrix, proj.method = "umap", dist.metric = "cosine", n.cores = 3, k, k.range,
                              var.genes = rownames(data.matrix), loss = "mse", genes.embed, hide.factors = T, n_pull = 3,
                              alpha.exp = 1.25, # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
                              snn.exp = 1.0 # Lower this < 1.0 to move similar cells closer to each other
){
  print(paste(length(var.genes), "variable genes"))
  if (missing(k)) {
    if (missing(k.range)) k.range <- seq(2,20,2) ## Range of factors to iterate over
    k.res <- FindNumFactors(data.matrix[var.genes,], k.range = k.range, n.cores = n.cores, do.plot = F, loss = loss)
    print(paste(k.res$k, "factors")); k <- k.res$k;
  }

  if (k < 3) {
    warning("k must be an integer of 3 or higher")
    k <- 3
  }

  pca.res <- irlba::irlba(t(data.matrix[var.genes,]), nv = max(k,20), center = Matrix::rowMeans(data.matrix[var.genes,]))
  pc.scores <- t(pca.res$u); colnames(pc.scores) <- colnames(data.matrix);
  snn <- CalcSNN(pc.scores, k = 20, prune.SNN = 1/20)

  if (missing(genes.embed)) genes.embed <- NULL
  run_swne(data.matrix, var.genes, snn, k, alpha.exp, snn.exp, n_pull, proj.method, dist.metric, genes.embed,
           loss, n.cores, hide.factors)
}



#' @rdname RunSWNE
#' @method RunSWNE matrix
#' @export
RunSWNE.matrix <- function(data.matrix, proj.method = "umap", dist.metric = "cosine", n.cores = 3, k, k.range,
                           var.genes = rownames(data.matrix), loss = "mse", genes.embed, hide.factors = T, n_pull = 3,
                           alpha.exp = 1.25, # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
                           snn.exp = 1.0 # Lower this < 1.0 to move similar cells closer to each other
){
  data.matrix <- as(data.matrix, "dgCMatrix")
  RunSWNE.dgCMatrix(data.matrix, proj.method = proj.method, dist.metric = dist.metric, n.cores = n.cores, k = k,
                    k.range = k.range, var.genes = var.genes, loss = loss, genes.embed = genes.embed,
                    hide.factors = hide.factors, n_pull = n_pull, alpha.exp = alpha.exp, snn.exp = snn.exp)
}



#' @rdname RunSWNE
#' @method RunSWNE dgTMatrix
#' @export
RunSWNE.dgTMatrix <- function(data.matrix, proj.method = "umap", dist.metric = "cosine", n.cores = 3, k, k.range,
                              var.genes = rownames(data.matrix), loss = "mse", genes.embed, hide.factors = T, n_pull = 3,
                              alpha.exp = 1.25, # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
                              snn.exp = 1.0 # Lower this < 1.0 to move similar cells closer to each other
){
  data.matrix <- as(data.matrix, "dgCMatrix")
  RunSWNE.dgCMatrix(data.matrix, proj.method = proj.method, dist.metric = dist.metric, n.cores = n.cores, k = k,
                    k.range = k.range, var.genes = var.genes, loss = loss, genes.embed = genes.embed,
                    hide.factors = hide.factors, n_pull = n_pull, alpha.exp = alpha.exp, snn.exp = snn.exp)
}



## Helper function for running SWNE
run_swne <- function(norm_counts, var.genes, snn, k, alpha.exp, snn.exp, n_pull, proj.method, dist.metric,
                     genes.embed, loss, n.cores, hide.factors) {
  nmf.res <- RunNMF(norm_counts[var.genes,], k = k, init = "ica", n.cores = n.cores, loss = loss)
  nmf.scores <- nmf.res$H
  swne_embedding <- EmbedSWNE(nmf.scores, snn, alpha.exp = alpha.exp, snn.exp = snn.exp,
                              n_pull = n_pull, proj.method = proj.method,
                              dist.use = dist.metric)

  if (!is.null(genes.embed)) {
    genes.embed <- intersect(genes.embed, rownames(norm_counts))
    nmf.loadings <- ProjectFeatures(norm_counts, nmf.scores, loss = loss, n.cores = n.cores)
    swne_embedding <- EmbedFeatures(swne_embedding, nmf.loadings, genes.embed, n_pull = n_pull)
  }
  if (hide.factors) swne_embedding$H.coords$name <- ""

  return(swne_embedding)
}
