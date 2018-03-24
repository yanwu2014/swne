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
#' @import FNN
#' @importFrom Matrix sparseMatrix
#'
#' @return Returns similarity matrix in sparse matrix format
#'
#' @import Seurat
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

  w <- compute_snn(train.cell.names = colnames(data.use), k = k, train.nn.large = nn.large,
                    train.nn.ranked = nn.ranked, prune.SNN = prune.SNN,
                    print.output = print.output)
  Matrix::diag(w) <- 1
  # w <- Seurat:::CalcSNNSparse(cell.names = colnames(data.use), k = k, nn.large = nn.large,
  #                             nn.ranked = nn.ranked, prune.SNN = prune.SNN,
  #                             print.output = print.output)
  return(w)
}


#' Calculates shared nearest neighbors between samples (columns) in a test matrix and samples in a training matrix
#' Adapted from Seurat
#'
#' @param test.matrix Test data
#' @param train.matrix Training data
#' @param k Number of nearest neighbors
#' @param k.scale k*k.scale is the number of nearest neighbors to calculate shared nearest neighbors for
#' @param prune.SNN Minimum fraction of shared nearest neighbors
#' @param print.output Prints progress bar if true
#'
#' @return A n.test x n.train matrix with the shared nearest neighbors between test and training data
#'
#' @export
#'
ProjectSNN <- function(test.matrix, train.matrix, k = 20, k.scale = 10, prune.SNN = 1/15, print.output = T) {
  n.train.cells <- ncol(train.matrix)
  n.test.cells <- ncol(test.matrix)
  stopifnot(k*k.scale < n.train.cells - 1)

  train.knn <- FNN::get.knn(data = t(train.matrix), k = k.scale * k)
  train.nn.ranked <- cbind(1:n.train.cells, train.knn$nn.index[, 1:(k - 1)])
  train.nn.large <- train.knn$nn.index

  test.knn <- FNN::get.knnx(data = t(train.matrix), query = t(test.matrix), k = k * k.scale)
  test.nn.ranked <- cbind(1:n.test.cells, test.knn$nn.index[, 1:(k - 1)])
  test.nn.large <- test.knn$nn.index

  w <- compute_snn(colnames(train.matrix), k, train.nn.large, train.nn.ranked,
                    test.cell.names = colnames(test.matrix),
                    test.nn.large = test.nn.large, test.nn.ranked = test.nn.ranked,
                    prune.SNN = prune.SNN, print.output = print.output)
  return(w)
}


#' Helper function for calculating a SNN graph
#' Adapted from Seurat
compute_snn <- function(train.cell.names, k, train.nn.large, train.nn.ranked, test.cell.names = NULL,
                        test.nn.large = NULL, test.nn.ranked = NULL, prune.SNN = 1/15, print.output = T) {

  if(any(is.null(c(test.cell.names, test.nn.large, test.nn.ranked)))) {
    print("Computing single dataset SNN")
    test.cell.names <- train.cell.names
    test.nn.large <- train.nn.large
    test.nn.ranked <- train.nn.ranked
  } else {
    print("Computing test data SNN projected onto training data")
  }

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
}
