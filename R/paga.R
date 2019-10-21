## Functions for plotting and using PAGA graphs to prune SNNs


#' Make PAGA graph. Adapted from Monocle3
#'
#' @param knn kNN in the form of a sparse matrix
#' @param clusters Cell clusters
#' @param qval.cutoff Q-value cutoff for cluster to cluster interactions
#'
#' @return List representing a PAGA graph
#'
#' @export
#'
BuildPAGA <- function(knn, clusters = NULL, qval.cutoff = 0.05) {
  if (!requireNamespace("monocle3", quietly = T)) {
    stop("Monocle3 is needed for this function to work. Please install it",
         call. = F)
  }

  if (!requireNamespace("igraph", quietly = T)) {
    stop("igraph is needed for this function to work. Please install it",
         call. = F)
  }

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
