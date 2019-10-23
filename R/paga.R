## Functions for plotting and using PAGA graphs to prune SNNs


#' Create a graph of cosine similarities between clusters pruned using PAGA
#'
#' @param partitions List representing a PAGA graph (output of BuildPAGA)
#' @param embeddings Dimensional reduction used to compute cosine similarities (cells x dimensions)
#' @param qval.cutoff Q-value cutoff for cluster to cluster interactions
#' @param min.cluster.sim Minimum similarity needed to draw an edge between two clusters
#' @param seed Random seed for force-directed layout reproducibility
#'
#' @return igraph object where each vertex is a cluster and edges represent PAGA-pruned similarities between clusters
#'
#' @import igraph
#' @export
#'
PartitionSimilarityGraph <- function(partitions, embedding, qval.cutoff = 1e-3,
                                     min.cluster.sim = 0.05, seed = 42) {
  stopifnot(all(rownames(embedding) %in% names(partitions$clusters)))
  embedding <- embedding[names(partitions$clusters),]
  mean.embedding <- apply(embedding, 2, function(x) tapply(x, partitions$clusters, mean))

  cluster.sim <- as.matrix(proxy::simil(mean.embedding, method = "cosine", by_rows = T))
  cluster.sim <- cluster.sim[rownames(partitions$cluster_mat), colnames(partitions$cluster_mat)]
  cluster.sim[partitions$cluster_mat > qval.cutoff] <- 0
  cluster.sim[cluster.sim < min.cluster.sim] <- 0

  cluster.graph <- graph_from_adjacency_matrix(cluster.sim, mode = "undirected", weighted = T)

  set.seed(seed)
  cluster.graph.layout <- layout_with_fr(cluster.graph, weights = E(cluster.graph)$weight)
  cluster.graph <- set_graph_attr(cluster.graph, "layout", cluster.graph.layout)

  return(cluster.graph)
}


#' Make PAGA graph. Adapted from Monocle3
#'
#' @param knn kNN in the form of a sparse matrix
#' @param clusters Cell clusters
#' @param qval.cutoff Q-value cutoff for cluster to cluster interactions
#'
#' @return List representing a PAGA graph
#'
#' @import igraph
#' @export
#'
BuildPAGA <- function(knn, clusters = NULL, qval.cutoff = 0.05) {
  knn <- as(knn, "dgCMatrix")
  knn_graph <- igraph::graph_from_adjacency_matrix(knn, mode = "undirected")

  knn_clusters <- igraph::cluster_louvain(knn_graph)
  if (is.null(clusters)) {
    knn_clusters$membership <- factor(knn_clusters$membership)
    table(knn_clusters$membership)
  } else {
    knn_clusters$membership <- factor(clusters)
  }

  names(knn_clusters$membership) <- names(igraph::V(knn_graph))
  partitions <- monocle3:::compute_partitions(knn_graph, knn_clusters,
                                              qval_thresh = qval.cutoff,
                                              verbose = F)

  partitions$clusters <- knn_clusters$membership
  rownames(partitions$cluster_mat) <- colnames(partitions$cluster_mat) <- levels(knn_clusters$membership)
  rownames(partitions$num_links) <- colnames(partitions$num_links) <- levels(knn_clusters$membership)

  return(partitions)
}



#' Prune SNN with PAGA graph. Adapted from Monocle3
#'
#' @param snn Shared Nearest Neighbors (SNN) graph in the form of a sparse matrix
#' @param knn kNN in the form of a sparse matrix
#' @param clusters Cell clusters
#' @param qval.cutoff Q-value cutoff for significant edges
#'
#' @return SNN graph with spurious edges pruned
#'
#' @export
#'
PruneSNN <- function(snn, knn, clusters = NULL, qval.cutoff = 0.05) {
  partitions <- BuildPAGA(knn, clusters = clusters, qval.cutoff = qval.cutoff)

  cluster_mat <- partitions$cluster_mat
  clusters <- partitions$clusters

  cl.cells.list <- lapply(levels(clusters), function(cl) {
    names(clusters[clusters == cl])
  })
  names(cl.cells.list) <- levels(clusters)

  snn <- as.matrix(snn)
  cl.combos <- combn(levels(clusters), 2)
  for(i in 1:ncol(cl.combos)) {
    cl1 <- cl.combos[1,i]
    cl2 <- cl.combos[2,i]
    if (cluster_mat[cl1, cl2] > qval.cutoff) {
      cl1.ix <- cl.cells.list[[cl1]]
      cl2.ix <- cl.cells.list[[cl2]]
      snn[cl1.ix, cl2.ix] <- snn[cl2.ix, cl1.ix] <- 0
    }
  }
  diag(snn) <- 1
  snn <- as(snn, "dgCMatrix"); invisible(gc());

  return(snn)
}



#' Compute significant links between specified clusters
#' Adapted from monocle3 (https://github.com/cole-trapnell-lab/monocle3)
#'
compute_partitions <- function(g,
                               optim_res,
                               qval_thresh=0.05,
                               verbose = FALSE){
  cell_membership <- as.factor(igraph::membership(optim_res))
  membership_matrix <- Matrix::sparse.model.matrix( ~ cell_membership + 0)
  num_links <- Matrix::t(membership_matrix) %*%
    igraph::as_adjacency_matrix(g) %*% membership_matrix
  diag(num_links) <- 0
  louvain_modules <- levels(cell_membership)

  edges_per_module <- Matrix::rowSums(num_links)
  total_edges <- sum(num_links)

  theta <- (as.matrix(edges_per_module) / total_edges) %*%
    Matrix::t(edges_per_module / total_edges)
  var_null_num_links <- theta * (1 - theta) / total_edges
  num_links_ij <- num_links / total_edges - theta
  cluster_mat <- pnorm_over_mat(as.matrix(num_links_ij), var_null_num_links)

  num_links <- num_links_ij / total_edges

  cluster_mat <- matrix(stats::p.adjust(cluster_mat),
                        nrow=length(louvain_modules),
                        ncol=length(louvain_modules))

  sig_links <- as.matrix(num_links)
  sig_links[cluster_mat > qval_thresh] = 0
  diag(sig_links) <- 0

  cluster_g <- igraph::graph_from_adjacency_matrix(sig_links, weighted = T,
                                                   mode = 'undirected')

  list(cluster_g = cluster_g, num_links = num_links, cluster_mat = cluster_mat)
}
