#### SNN functions ####

#' SNN Graph Construction. Adapted from Seurat.
#'
#' @param data.use Features x samples matrix to use to build the SNN
#' @param k Defines k for the k-nearest neighbor algorithm
#' @param k.scale Granularity option for k.param
#' @param prune.SNN Sets the cutoff for acceptable Jaccard distances when
#'                  computing the neighborhood overlap for the SNN construction.
#' @param print.output Whether or not to print output to the console
#'
#' @return Returns similarity matrix in sparse matrix format
#'
#' @importFrom FNN get.knn
#' @importFrom Matrix sparseMatrix
#' @export
#'
CalcSNN <- function(data.use, k = 10, k.scale = 10, prune.SNN = 1/15, print.output = T) {
  n.cells <- ncol(data.use)
  if (n.cells < k) {
    stop("k cannot be greater than the number of samples")
  }

  ## find the k-nearest neighbors for each single cell
  my.knn <- FNN::get.knn(t(as.matrix(data.use)), k = min(k.scale * k, n.cells - 1))
  nn.ranked <- cbind(1:n.cells, my.knn$nn.index[, 1:(k - 1)])
  nn.large <- my.knn$nn.index

  w <- ComputeSNN(nn.ranked, prune.SNN)
  colnames(w) <- rownames(w) <- colnames(data.use)

  Matrix::diag(w) <- 1
  return(w)
}


#' Calculates shared nearest neighbors between samples (columns) in a test matrix and samples in a training matrix
#' Adapted from Seurat
#'
#' @param test.data Test data
#' @param train.data Training data
#' @param pcs.use If not null, run PCA on training data and project test data before computing SNN
#' @param features.use Subset of features to use (default is all features)
#' @param k Number of nearest neighbors
#' @param k.scale k*k.scale is the number of nearest neighbors to calculate shared nearest neighbors for
#' @param prune.SNN Minimum fraction of shared nearest neighbors
#' @param print.output Prints progress bar if true
#'
#' @return A n.test x n.train matrix with the shared nearest neighbors between test and training data
#'
#' @importFrom FNN get.knn
#' @importFrom Matrix sparseMatrix
#' @import irlba
#' @export
#'
ProjectSNN <- function(test.data, train.data, n.pcs = NULL, features.use = NULL, k = 30, k.scale = 10,
                       prune.SNN = 1/15, print.output = T) {
  n.train.cells <- ncol(train.data)
  n.test.cells <- ncol(test.data)
  stopifnot(k*k.scale < n.train.cells - 1)

  if (is.null(features.use)) {
    features.use <- intersect(rownames(test.data), rownames(train.data))
  }

  if (!is.null(n.pcs)) {
    print("Running PCA and test data projection")
    train.cm <- Matrix::rowMeans(train.data[features.use,])
    train.pca <- irlba::irlba(t(train.data[features.use,] - train.cm), nv = n.pcs)
    train.loadings <- train.pca$v; rownames(train.loadings) <- features.use;

    test.cm <- Matrix::rowMeans(test.data[features.use,])
    test.data.use <- t(t(test.data[features.use,] - test.cm) %*% train.loadings)
    train.data.use <- t(train.pca$u); colnames(train.data.use) <- colnames(train.data);
  } else {
    test.data.use <- test.data[features.use,]
    train.data.use <- train.data[features.use,]
  }


  train.knn <- FNN::get.knn(data = t(train.data.use), k = k.scale * k)
  train.nn.ranked <- cbind(1:n.train.cells, train.knn$nn.index[, 1:(k - 1)])

  test.knn <- FNN::get.knnx(data = t(train.data.use), query = t(test.data.use), k = k * k.scale)
  test.nn.ranked <- cbind(1:n.test.cells, test.knn$nn.index[, 1:(k - 1)])
  test.nn.large <- test.knn$nn.index

  w <- compute_projected_snn(colnames(train.data.use), k, train.nn.ranked, colnames(test.data),
                             test.nn.large, test.nn.ranked, prune.SNN, print.output)
  return(w)
}


#' Helper function for calculating a SNN graph
#' Adapted from Seurat
#'
#' @import compiler
#'
compute_projected_snn <- function(train.cell.names, k, train.nn.ranked, test.cell.names,
                                  test.nn.large, test.nn.ranked, prune.SNN = 1/15,
                                  print.output = T) {
  print("Computing SNN of test data projected onto training data")

  n.train.cells <- length(train.cell.names)
  n.test.cells <- length(test.cell.names)

  counter <- 1
  idx1 <- vector(mode = "integer", length = n.test.cells * n.train.cells / k)
  idx2 <- vector(mode = "integer", length = n.test.cells * n.train.cells / k)
  edge.weight <- vector(mode = "double", length = n.test.cells * n.train.cells / k)
  id <- 1

  if (print.output) {
    print("Constructing SNN")
    pb <- txtProgressBar(min = 0, max = n.test.cells, style = 3)
  }

  for (i in 1:n.test.cells) {
    for (j in 1:ncol(test.nn.large)) {
      s <- intersect(test.nn.ranked[i, ], train.nn.ranked[test.nn.large[i, j], ])
      u <- union(test.nn.ranked[i, ], train.nn.ranked[test.nn.large[i, j], ])
      e <- length(s) / length(u)
      if (e > prune.SNN) {
        idx1[id] <- i
        idx2[id] <- test.nn.large[i, j]
        edge.weight[id] <- e
        id <- id + 1
      }
    }
    if (print.output) {
      setTxtProgressBar(pb = pb, value = i)
    }
  }
  if (print.output) {
    close(con = pb)
  }

  idx1 <- idx1[! is.na(idx1) & idx1 != 0]
  idx2 <- idx2[! is.na(idx2) & idx2 != 0]
  edge.weight <- edge.weight[! is.na(edge.weight) & edge.weight != 0]
  w <- sparseMatrix(
    i = idx1,
    j = idx2,
    x = edge.weight,
    dims = c(n.test.cells, n.train.cells)
  )

  rownames(w) <- test.cell.names
  colnames(w) <- train.cell.names
  return(w)
}; compute_projected_snn <- compiler::cmpfun(compute_projected_snn)

