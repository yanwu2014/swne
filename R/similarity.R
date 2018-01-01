## SNN functions

.calc_snn_project <- function(train.cell.names, test.cell.names, k, test.nn.large, test.nn.ranked, train.nn.ranked,
                              prune.SNN, print.output = T) {
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
    for (j in 1:ncol(x = test.nn.large)) {
      s <- intersect(x = test.nn.ranked[i, ], y = train.nn.ranked[test.nn.large[i, j], ])
      u <- union(test.nn.ranked[i, ], train.nn.ranked[test.nn.large[i, j], ])
      e <- length(x = s) / length(x = u)
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

  idx1 <- idx1[! is.na(x = idx1) & idx1 != 0]
  idx2 <- idx2[! is.na(x = idx2) & idx2 != 0]
  edge.weight <- edge.weight[! is.na(x = edge.weight) & edge.weight != 0]
  w <- sparseMatrix(
    i = idx1,
    j = idx2,
    x = edge.weight,
    dims = c(n.test.cells, n.train.cells)
  )

  rownames(x = w) <- test.cell.names
  colnames(x = w) <- train.cell.names
  return(w)
}



ProjectSNN <- function(test.matrix, train.matrix, k = 20, k.scale = 10, prune.SNN = 1/15, print.output = T) {
  n.train.cells <- ncol(train.matrix)
  n.test.cells <- ncol(test.matrix)
  stopifnot(k*k.scale < n.train.cells - 1)

  train.knn <- FNN::get.knn(data = t(train.matrix), k = k.scale * k)
  train.nn.ranked <- cbind(1:n.train.cells, train.knn$nn.index[, 1:(k - 1)])

  test.knn <- FNN::get.knnx(data = t(train.matrix), query = t(test.matrix), k = k * k.scale)
  test.nn.ranked <- cbind(1:n.test.cells, test.knn$nn.index[, 1:(k - 1)])
  test.nn.large <- test.knn$nn.index

  test.snn <- .calc_snn_project(colnames(train.matrix), colnames(test.matrix), k, test.nn.large, test.nn.ranked, train.nn.ranked,
                                prune.SNN, print.output = print.output)
  return(test.snn)
}
