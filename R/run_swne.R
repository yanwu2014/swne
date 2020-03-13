#' Wrapper for running SWNE analysis
#'
#' @param object A Seurat, Pagoda2 or cisTopicObject object with normalized data
#' @param reduction.use Which dimensional reduction (e.g. PCA, ICA) to use for the tSNE. Default is PCA.
#' @param cells.use Which cells to analyze (default, all cells)
#' @param dims.use Which dimensions to use as input features
#' @param genes.use If set, run the SWNE on this subset of genes (instead of running on a set of reduced dimensions). Not set (NULL) by default
#' @param data.matrix a data matrix (genes x cells) which has been pre-normalized
#' @param batch Vector of batch IDs to regress away
#' @param proj.method Method to use to project factors in 2D. Either "sammon" or "umap"
#' @param dist.use Similarity function to use for calculating factor positions (passed to EmbedSWNE).
#'                 Options include pearson (correlation), IC (mutual information), cosine, euclidean.
#' @param distance.matrix If set, runs tSNE on the given distance matrix instead of data matrix (experimental)
#' @param n.cores Number of cores to use (passed to FindNumFactors)
#' @param k Number of NMF factors (passed to RunNMF). If none given, will be derived from FindNumFactors.
#' @param k.range Range of factors for FindNumFactors to iterate over if k is not given
#' @param var.genes vector to specify variable genes. Will infer from Seurat or use full dataset if not given.
#' @param loss loss function to use (passed to RunNMF)
#' @param alpha.exp Increasing alpha.exp increases how much the NMF factors "pull" the samples (passed to EmbedSWNE)
#' @param snn.exp Decreasing snn.exp increases the effect of the similarity matrix on the embedding (passed to EmbedSWNE)
#' @param snn.k Changes the number of nearest neighbors used to build SNN (passed to CalcSNN)
#' @param prune.SNN The minimum fraction of shared nearest neighbors (smaller values are set to zero)
#' @param use.paga.pruning Use PAGA graphs to prune
#' @param sample.groups Clusters to use for PAGA (default is to do a de-novo clustering)
#' @param paga.qval.cutoff q-value cutoff for significant shared edges between clusters
#' @param n.var.genes Number of variable genes to use
#' @param n_pull Maximum number of factors "pulling" on each sample
#' @param ica.fast Whether to run SVD before ICA initialization
#' @param genes.embed Genes to add to the SWNE embedding
#' @param hide.factors Hide factors when plotting SWNE embedding
#' @param reduction.name dimensional reduction name, specifies the position in the object$dr list. swne by default
#' @param reduction.key dimensional reduction key, specifies the string before the number for the dimension names. SWNE_ by default
#' @param return.format format to return ("seurat" object or raw "embedding")
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
#' @method RunSWNE cisTopic
#' @export
RunSWNE.cisTopic <- function(cisTopicObject, proj.method = "sammon", cells.use = NULL,
                             dist.metric = "cosine", n.cores = 8, hide.factors = T, n_pull = 3,
                             alpha.exp = 1.25, # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
                             snn.exp = 1.0, # Lower this < 1.0 to move similar cells closer to each other
                             snn.k = 20,
                             prune.SNN = 1/15,
                             use.paga.pruning = T,
                             sample.groups = NULL,
                             paga.qval.cutoff = 1e-3) {
  if (!requireNamespace("cisTopic", quietly = T)) {
    stop("cisTopic needed for this function to work. Please install it.",
         call. = F)
  }

  cisTopicObject <- getRegionsScores(cisTopicObject)
  cisTopicObject <- cisTopic::runPCA(cisTopicObject, target = "cell", method = "Probability")

  pc.emb <- t(cisTopicObject@dr$cell[["PCA"]]$ind.coord)
  topic.emb <- modelMatSelection(cisTopicObject, target = "cell", method = "Probability")

  snn <- CalcSNN(pc.emb, k = snn.k, prune.SNN = prune.SNN)
  if (use.paga.pruning) {
    knn <- CalcKNN(pc.emb, k = snn.k)
    snn <- PruneSNN(snn, knn, clusters = sample.groups, qval.cutoff = paga.qval.cutoff)
  }
  swne.emb <- EmbedSWNE(topic.emb, snn, alpha.exp = alpha.exp, snn.exp = snn.exp, n_pull = n_pull)

  if (hide.factors) {
    swne.emb$H.coords$name <- ""
  }

  return(swne.emb)
}



#' @rdname RunSWNE
#' @method RunSWNE Seurat
#' @export
RunSWNE.Seurat <- function(object, proj.method = "sammon", reduction.use = "pca", cells.use = NULL, dims.use = NULL, genes.use = NULL,
                           dist.metric = "cosine", distance.matrix = NULL, n.cores = 8, k, k.range, var.genes,
                           loss = "mse", genes.embed, hide.factors = T, n_pull = 3, ica.fast = T,
                           alpha.exp = 1.25, # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
                           snn.exp = 1.0, # Lower this < 1.0 to move similar cells closer to each other
                           snn.k = 10,
                           use.paga.pruning = T,
                           sample.groups = NULL,
                           paga.qval.cutoff = 1e-3,
                           reduction.name = "swne", reduction.key = "SWNE_", return.format = "embedding", ...
){
  if (!requireNamespace("Seurat", quietly = T)) {
    stop("Seurat is needed for this function to work. Please install it",
         call. = F)
  }

  if (is.null(dims.use)) {
    if (is.null(k) || missing(k)) {
      dims.use <- 1:20
    } else {
      dims.use <- 1:k
    }
  }
  if (length(x = dims.use) < 3) {
    stop("SWNE needs at least 3 dimensions to generate a plot")
  }
  if (!is.null(x = distance.matrix)) {
    genes.use <- rownames(x = object)
  }

  if (DefaultAssay(object) == "integrated") {
    object_norm <- as.matrix(GetAssayData(object, assay = "integrated"))
    object_norm <- t(apply(object_norm, 1, function(x) (x - min(x))/(max(x) - min(x))))
    var.genes <- rownames(object_norm)
  } else {
    object_norm <- ExtractNormCounts(object, obj.type = "seurat", rescale = F, rescale.method = "log", batch = NULL)
  }

  if (missing(var.genes)) var.genes <- VariableFeatures(object)
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

  if (missing(genes.embed)) genes.embed <- NULL
  if (is.null(distance.matrix)) {
    if (DefaultAssay(object) == "integrated") {
      if(sum(dim(object@graphs$integrated_snn)) != 2*ncol(object)) {
        object <- RunPCA(object, pc.genes = var.genes, do.print = F, pcs.compute = min(k,20),
                         verbose = F)
        object <- FindNeighbors(object, k = snn.k, prune.SNN = 1/15)
      }
      snn <- as(object@graphs$integrated_snn, "dgCMatrix")
      if (use.paga.pruning) knn <- as(object@graphs$integrated_nn, "dgCMatrix")
    } else {
      if(sum(dim(object@graphs$RNA_snn)) != 2*ncol(object)) {
        object <- RunPCA(object, pc.genes = var.genes, do.print = F, pcs.compute = min(k,20),
                         verbose = F)
        object <- FindNeighbors(object, k = snn.k, prune.SNN = 1/15)
      }
      snn <- as(object@graphs$RNA_snn, "dgCMatrix")
      if (use.paga.pruning) knn <- as(object@graphs$RNA_nn, "dgCMatrix")
    }

    if (use.paga.pruning) snn <- PruneSNN(snn, knn, clusters = sample.groups, qval.cutoff = paga.qval.cutoff)
    swne_embedding <- run_swne(object_norm, var.genes, snn, k, alpha.exp, snn.exp, n_pull, proj.method, dist.metric, genes.embed,
                               loss, n.cores, hide.factors, ica.fast)
  }
  else {
    swne_embedding <- RunSWNE(as.matrix(distance.matrix), proj.method = proj.method, dist.metric = dist.metric, n.cores = n.cores, k = k,
                              k.range = k.range, var.genes = var.genes, loss = loss, genes.embed = genes.embed,
                              hide.factors = hide.factors, n_pull = n_pull, alpha.exp = alpha.exp, snn.exp = snn.exp)
  }

  if(return.format == "embedding"){
    return(swne_embedding)
  } else if(return.format == "seurat"){
    swne.emb.scores <- as.matrix(swne_embedding$sample.coords)
    colnames(swne.emb.scores) <- paste0(reduction.key, 1:ncol(swne.emb.scores))
    swne_dr <- CreateDimReducObject(embeddings = swne.emb.scores,
                                    key = reduction.key, assay = DefaultAssay(object),
                                    misc = list(genes.use = var.genes,
                                                cells.use = colnames(object)))
    object[[reduction.name]] <- swne_dr
    return(object)
  }
}



#' @rdname RunSWNE
#' @method RunSWNE Pagoda2
#' @export
#'
RunSWNE.Pagoda2 <- function(object, proj.method = "sammon", dist.metric = "cosine", n.cores = 8, k, k.range, var.genes,
                            loss = "mse", genes.embed, hide.factors = T, n_pull = 3, ica.fast = T, n.var.genes = 3000,
                            alpha.exp = 1.25, # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
                            snn.exp = 1.0, # Lower this < 1.0 to move similar cells closer to each other
                            snn.k = 20,
                            use.paga.pruning = T,
                            sample.groups = NULL,
                            paga.qval.cutoff = 1e-3
){
  if (!requireNamespace("pagoda2", quietly = T)) {
    stop("pagoda2 needed for this function to work. Please install it.",
         call. = F)
  }

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
  snn <- CalcSNN(pc.scores, k = snn.k, prune.SNN = 1/20)

  if (use.paga.pruning) {
    knn <- CalcKNN(pc.scores, k = snn.k)
    snn <- PruneSNN(snn, knn, clusters = sample.groups, qval.cutoff = paga.qval.cutoff)
  }

  if (missing(genes.embed)) genes.embed <- NULL
  run_swne(object_norm, var.genes, snn, k, alpha.exp, snn.exp, n_pull, proj.method, dist.metric, genes.embed,
           loss, n.cores, hide.factors, ica.fast)
}



#' @rdname RunSWNE
#' @method RunSWNE dgCMatrix
#' @export
RunSWNE.dgCMatrix <- function(data.matrix, proj.method = "sammon", dist.metric = "cosine", n.cores = 3, k, k.range,
                              var.genes, loss = "mse", genes.embed, hide.factors = T,
                              n_pull = 3, ica.fast = T,
                              alpha.exp = 1.25, # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
                              snn.exp = 1.0, # Lower this < 1.0 to move similar cells closer to each other
                              snn.k = 20,
                              use.paga.pruning = T,
                              sample.groups = NULL,
                              paga.qval.cutoff = 1e-3
){
  if (missing(var.genes)) {
    var.genes = rownames(data.matrix)
  }
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
  snn <- CalcSNN(pc.scores, k = snn.k, prune.SNN = 1/20)

  if (use.paga.pruning) {
    knn <- CalcKNN(pc.scores, k = snn.k)
    snn <- PruneSNN(snn, knn, clusters = sample.groups, qval.cutoff = paga.qval.cutoff)
  }

  if (missing(genes.embed)) genes.embed <- NULL
  run_swne(data.matrix, var.genes, snn, k, alpha.exp, snn.exp, n_pull, proj.method, dist.metric, genes.embed,
           loss, n.cores, hide.factors, ica.fast)
}



#' @rdname RunSWNE
#' @method RunSWNE matrix
#' @export
RunSWNE.matrix <- function(data.matrix, proj.method = "sammon", dist.metric = "cosine", n.cores = 8, k, k.range,
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
RunSWNE.dgTMatrix <- function(data.matrix, proj.method = "sammon", dist.metric = "cosine", n.cores = 3, k, k.range,
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
                     genes.embed, loss, n.cores, hide.factors, ica.fast) {
  nmf.res <- RunNMF(norm_counts[var.genes,], k = k, init = "ica", n.cores = n.cores, loss = loss,
                    ica.fast = ica.fast)
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
