## Geneset handling functions

#' Load genesets from gmt file. For more info on the .gmt file format, see [link]
#'
#' @param geneset.file Geneset filename
#'
#' @return List of genesets (character vectors)
#'
#' @export
#'
ReadGenesets <- function(genesets.file) {
  genesets <- as.list(readLines(genesets.file))
  genesets <- lapply(genesets, function (v) strsplit(v, '\t')[[1]])
  geneset.names <- unlist(lapply(genesets, function(x) x[[1]]))
  genesets <- lapply(genesets, function(v) v[3:length(v)])
  names(genesets) <- sapply(geneset.names, function(x) gsub(" ", "_", x))
  return(genesets)
}


## Helper function for filtering genesets
.clean_genesets <- function(genesets, min.size = 5, max.size = 500, annot = FALSE) {
  genesets <- as.list(genesets)
  size <- unlist(lapply(genesets, length))
  genesets <- genesets[size > min.size & size < max.size]
  return(genesets)
}


#' Filters genesets according to minimum and maximum set size
#'
#' @param gene.names Name of all genes to use
#' @param min.size Minimum geneset size
#' @param max.size Maximum geneset size
#'
#' @return Filtered list of genesets
#'
#' @export
#'
FilterGenesets <- function(genesets, gene.names, min.size = 5, max.size = 500) {
  genesets <- lapply(genesets, function (x) return(x[x %in% gene.names]))
  genesets <- .clean_genesets(genesets, min.size = min.size, max.size = max.size)
  return(genesets)
}


#' Writes genesets to gmt file
#'
#' @param genesets List of genesets
#' @param file.name Name of output file
#'
#' @return Writes to file
#'
#' @export
#'
WriteGenesets <- function(genesets, file.name) {
  genesets <- lapply(names(genesets), function(name) {x <- genesets[[name]]; x <- c(name,name,x); return(x);})
  n.cols <- 1.5*max(unlist(lapply(genesets, length)))
  empty <- lapply(genesets, write, file.name, append=TRUE, ncolumns = n.cols, sep = '\t')
}


## Functions for annotating NMF factors

## Creates a genes x genesets indicator matrix
.genesets_indicator <- function(genesets, inv = F, return.numeric = F) {
  genes <- unique(unlist(genesets, F, F))

  ind <- matrix(F, length(genes), length(genesets))
  rownames(ind) <- genes
  colnames(ind) <- names(genesets)

  for (i in 1:length(genesets)) {
    ind[genesets[[i]], i] <- T
  }

  if(inv) {
    ind <- !ind
  }

  if (return.numeric) {
    ind.numeric <- apply(ind, 2, function(x) as.numeric(x))
    rownames(ind.numeric) <- rownames(ind)
    return(ind.numeric)
  } else {
    return(ind)
  }
}


#' Projects dataset onto genesets using either NMF or a nonnegative linear model
#'
#' @param norm.counts Normalized input data
#' @param genesets List of genesets to use
#' @param method Method used to project data onto genesets. lm uses a simple linear model while mf
#'               will use a masked version of nmf.
#' @param loss Loss function for the nonnegative linear model
#'
#' @return Genesets x samples matrix
#'
#' @export
#'
ProjectGenesets <- function(norm.counts, genesets, method = "lm", loss = "mse") {
  if (!method %in% c("lm", "mf")) { stop("Invalid method") }

  if (method == "lm") {
    full.genesets.matrix <- .genesets_indicator(genesets, inv = F, return.numeric = T)
    genesets.H <- ProjectNMF(norm.counts[rownames(full.genesets.matrix),], full.genesets.matrix,
                                  loss = loss)
  } else if (method == "mf") {
    genesets.mask <- .genesets_indicator(genesets, inv = T)
    nmf.res <- nnmf(as.matrix(norm.counts[rownames(genesets.mask),]), k = ncol(genesets.mask), loss = loss,
                    mask = list(W = genesets.mask))
    colnames(nmf.res$W) <- rownames(nmf.res$H) <- colnames(genesets.mask)
    genesets.H <- nmf.res$H
  }

  genesets.H
}


#' Find the top gene or geneset markers for each NMF factor using pearson/spearman correlation,
#' or mutual information
#'
#' @param feature.mat Feature matrix (features x samples)
#' @param H Factor scores (factors x samples)
#' @param n.cores Number of cores to use
#' @param metric Association metric: pearson, spearman, or IC (information coefficient)
#'
#' @return Features x factors matrix of associations or correlations
#'
#' @export
#'
FactorAssociation <- function(feature.mat, nmf.scores, n.cores = 8, metric = "IC") {
  cl <- snow::makeCluster(n.cores, type = "SOCK")
  snow::clusterExport(cl, c("nmf.scores", "MutualInf"), envir = environment())
  if (metric == "IC") {
    source.log <- snow::parLapply(cl, 1:length(cl), function(i) library(MASS))
    assoc <- t(snow::parApply(cl, feature.mat, 1, function(v) apply(nmf.scores, 1, function(u) MutualInf(u, v))))
  } else if (metric == "pearson") {
    assoc <- t(snow::parApply(cl, feature.mat, 1, function(v) apply(nmf.scores, 1, function(u) cor(u, v))))
  } else if (metric == "spearman") {
    assoc <- t(snow::parApply(cl, feature.mat, 1, function(v) apply(nmf.scores, 1, function(u) cor(u, v, method = "spearman"))))
  } else {
    stop("Invalid correlation metric")
  }
  stopCluster(cl)

  rownames(assoc) <- rownames(feature.mat)
  colnames(assoc) <- rownames(nmf.scores)
  return(assoc)
}


#' Summarize the factors x features association matrix into a readable dataframe
#'
#' @param feature.factor.assoc (features x factors) association matrix from FactorAssociation
#' @param features.return Number of top features to return for each factor
#' @param features.use Only consider a subset of features. Default is NULL, which will use all features
#'
#' @return Dataframe summarizing top features for each factor
#'
#' @export
#'
SummarizeAssocFeatures <- function(feature.factor.assoc, features.return = 10, features.use = NULL) {
  if (!is.null(features.use)) {
    feature.factor.assoc <- feature.factor.assoc[features.use,]
  }

  factor.features.df <- do.call("rbind", lapply(1:ncol(feature.factor.assoc), function(i) {
    features.df <- data.frame(assoc_score = feature.factor.assoc[,i])
    features.df$feature <- rownames(feature.factor.assoc)
    features.df$factor <- colnames(feature.factor.assoc)[[i]]
    features.df <- features.df[order(features.df$assoc_score, decreasing = T),]
    head(features.df, n = features.return)
  }))

  rownames(factor.features.df) <- NULL
  return(factor.features.df)
}


#' Runs GSEA on the gene association coefficients for each NMF
#'
#' @param gene.factor.assoc Matrix with the gene associations for each nmf
#' @param genesets Genesets to use (as a named list)
#' @param power GSEA coefficient power
#' @param n.rand Number of permutations to use when calculating significance
#' @param n.cores Number of threads to use
#'
#' @return Dataframe summarizing the gsea results
#'
#' @import liger
#'
#' @export
#'
RunGSEA <- function(gene.factor.assoc, genesets, power = 1, n.rand = 1000, n.cores = 1) {
  nmfs <- colnames(gene.factor.assoc)
  gsea.list <- lapply(nmfs, function(nf) {
    gene.assocs <- gene.factor.assoc[,nf]
    gene.assocs <- sort(gene.assocs, decreasing = T)
    df <- liger::bulk.gsea(gene.assocs, genesets, power = power, n.rand = n.rand, mc.cores = n.cores)
    df$geneset <- rownames(df); df$nmf <- nf;
    df
  })

  gsea.df <- do.call("rbind", gsea.list); rownames(gsea.df) <- NULL;
  gsea.df
}


#' Compute Information Coefficient [IC]
#' Pablo Tamayo Dec 30, 2015
#'
#' @param x Input vector x
#' @param y Input vector y
#' @param n.grid Gridsize for calculating IC
#'
#' @return Mutual information between x and y
#'
#' @export
#'
MutualInf <-  function(x, y, n.grid = 25) {
  x.set <- !is.na(x)
  y.set <- !is.na(y)
  overlap <- x.set & y.set

  x <- x[overlap] +  0.000000001*runif(length(overlap))
  y <- y[overlap] +  0.000000001*runif(length(overlap))

  if (length(x) > 2) {
    delta = c(MASS::bcv(x), MASS::bcv(y))
    rho <- cor(x, y)
    rho2 <- abs(rho)
    delta <- delta*(1 + (-0.75)*rho2)
    kde2d.xy <- MASS::kde2d(x, y, n = n.grid, h = delta)
    FXY <- kde2d.xy$z + .Machine$double.eps
    dx <- kde2d.xy$x[2] - kde2d.xy$x[1]
    dy <- kde2d.xy$y[2] - kde2d.xy$y[1]
    PXY <- FXY/(sum(FXY)*dx*dy)
    PX <- rowSums(PXY)*dy
    PY <- colSums(PXY)*dx
    HXY <- -sum(PXY * log(PXY))*dx*dy
    HX <- -sum(PX * log(PX))*dx
    HY <- -sum(PY * log(PY))*dy
    PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
    PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
    MI <- sum(PXY * log(PXY/(PX*PY)))*dx*dy
    IC <- sign(rho) * sqrt(1 - exp(- 2 * MI))
    if (is.na(IC)) IC <- 0
  } else {
    IC <- 0
  }
  return(IC)
}


## Genotype/group handling functions

#' Convert list of sample groups into a flat vector.
#' Removes all samples belonging to multiple groups
#'
#' @param groups.list List of sample groups
#'
#' @return A character vector of the group each sample belongs to. The vector names are the sample names.
#'
#' @export
#'
FlattenGroups <- function(groups.list) {
  cell.names <- unlist(groups.list, F, F)
  if (length(cell.names) > length(unique(cell.names))) {
    print("Warning: removing all samples belonging to multiple groups")
    groups.list <- single.groups(groups.list)
    cell.names <- unlist(groups.list, F, F)
  }

  groups <- rep(NA, length(cell.names))
  names(groups) <- cell.names
  for (g in names(groups.list)) {
    groups[groups.list[[g]]] <- g
  }

  if (any(is.na(groups))) { stop("Some unassigned cells"); }
  return(groups)
}


#' Converts a flat sample groups character vector into a list format
#'
#' @param groups Character vector of sample groups
#'
#' @return List of groups: each list element contains the samples for that group
#'
#' @export
#'
UnflattenGroups <- function(groups) {
  groups.list <- c()
  unique.groups <- unique(groups)
  for (g in unique.groups) {
    g.cells <- names(groups[groups == g])
    groups.list[[g]] <- g.cells
  }
  return(groups.list)
}


#' Read in sample groups from a csv file format
#'
#' @param groups.file Name of groups file
#' @param sep.char Delimiter
#'
#' @return List of groups
#'
#' @export
#'
ReadGroups <- function(groups.file, sep.char = ",") {
  group.data <- read.table(groups.file, sep = sep.char, header = F, stringsAsFactors = F)
  groups <- group.data[[1]]
  full.groups.list <- lapply(rownames(group.data), function(i) sapply(strsplit(group.data[i,2], split = ',')[[1]], trimws))
  full.groups.list <- lapply(full.groups.list, function(x) make.names(x))
  names(full.groups.list) <- groups
  return(full.groups.list)
}


#' Write sample groups from list to csv format
#'
#' @param groups.list List of groups
#' @param out.file Name of output file
#'
#' @export
#'
WriteGroups <- function(groups.list, out.file) {
  group.data <- sapply(groups.list, function(x) paste('\"', paste(x, collapse = ", "), '\"', sep = ""))
  group.data <- sapply(names(group.data), function(x) paste(x, group.data[[x]], sep = ","))
  fileConn = file(out.file)
  writeLines(group.data, fileConn)
  close(fileConn)
}
