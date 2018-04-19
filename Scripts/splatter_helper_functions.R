library(FNN)

## Calculate pairwise distances between centroids
CalcPairwiseDist <- function(data.use, clusters, dist.method = "euclidean", use.median = F) {
  if (use.median) {
    data.centroids <- t(apply(data.use, 1, function(x) tapply(x, clusters, median)))
  } else {
    data.centroids <- t(apply(data.use, 1, function(x) tapply(x, clusters, mean)))
  }
  return(proxy::dist(data.centroids, method = dist.method, by_rows = F))
}; CalcPairwiseDist <- compiler::cmpfun(CalcPairwiseDist);


## Function for identifying cells in the same path and time step
GetPathStep <- function(metadata, step.size, make.factor = T) {
  path.step <- as.character(metadata$Group); names(path.step) <- rownames(metadata);
  for (path in levels(factor(metadata$Group))) {
    steps <- sort(unique(subset(metadata, Group == path)$Step))
    step.range <- seq(min(steps), max(steps), step.size)
    for(i in step.range) {
      cells.step <- rownames(subset(metadata, Group == path & Step %in% seq(i, i + step.size - 1, 1)))
      path.step[cells.step] <- paste(path, i, sep = ".")
    }
  }
  if (make.factor) {
    path.step <- factor(path.step)
  }
  path.step
}


## Calculate approximate kNN for an embedding
ComputeKNN <- function(emb, k) {
  knn.idx <- knn.index(t(emb), k = k)
  knn.matrix <- matrix(0, ncol(emb), ncol(emb))
  for (i in 1:nrow(knn.idx)) {
    knn.matrix[knn.idx[i,],i] <- 1
    knn.matrix[i, knn.idx[i,]] <- 1
  }
  rownames(knn.matrix) <- colnames(knn.matrix) <- colnames(emb)
  as(knn.matrix, "dgCMatrix")
}; ComputeKNN <- compiler::cmpfun(ComputeKNN);


## Calculate Jaccard similarities
CalcJaccard <- function(x,y) {
  a <- sum(x)
  b <- sum(y)
  c <- sum(x == 1 & y == 1)
  c/(a + b - c)
}; CalcJaccard <- compiler::cmpfun(CalcJaccard);


## Make barplot from vector
ggBarplot <- function(x, y.lim = NULL, fill.color = "tomato") {
  require(ggplot2)
  barplot.df <- data.frame(Y = x, X = factor(names(x), levels = names(x)))
  
  ggobj <- ggplot(data = barplot.df, aes(x = X, y = Y)) + 
    geom_bar(position = "dodge", stat = "identity", width = 0.8, fill = fill.color, color = "black") + 
    theme_classic() + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
          axis.text.x = element_text(hjust = 1, size = 14, angle = 90, color = "black"), 
          axis.text.y = element_text(size = 12, color = "black"))
  if(!is.null(y.lim)) {
    ggobj <- ggobj + coord_cartesian(ylim = y.lim)
  }
  ggobj
}
