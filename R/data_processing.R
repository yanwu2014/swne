#### Matrix processing functions

#' Extract a field from a delimited string
#'
#' @param string Input string
#' @param field Field to extract
#' @param delim Character delimiter
#'
#' @export
#'
ExtractField <- function (string, field = 1, delim = "_") {
  fields <- as.numeric(unlist(strsplit(x = as.character(x = field), split = ",")))
  if (length(fields) == 1) {
    return(strsplit(string, split = delim)[[1]][field])
  }
  return(paste(strsplit(string, split = delim)[[1]][fields], collapse = delim))
}


#' Reads data from a sparse or dense matrix
#' Adapted from Seurat
#'
#' @param matrix.dir The input matrix or sparse matrix directory
#' @param make.sparse Force output to be a sparse dgCMatrix
#'
#' @return Returns either a dense or sparse matrix
#'
#' @import Matrix
#' @export
#'
ReadData <- function(matrix.dir, make.sparse = T) {
  if (dir.exists(matrix.dir)) {
    if (!grepl("\\/$", matrix.dir)) { matrix.dir <- paste(matrix.dir, "/", sep = "") }

    barcode.loc <- paste0(matrix.dir, "barcodes.tsv")
    gene.loc <- paste0(matrix.dir, "genes.tsv")
    matrix.loc <- paste0(matrix.dir, "matrix.mtx")

    if (!file.exists(barcode.loc)) { stop("Barcode file missing") }
    if (!file.exists(gene.loc)) { stop("Gene name file missing") }
    if (!file.exists(matrix.loc)) { stop("Expression matrix file missing") }

    counts <- Matrix::readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    gene.names <- readLines(gene.loc)
    if (all(grepl(pattern = "\\-1$", cell.names))) {
      cell.names <- as.vector(as.character(sapply(cell.names, ExtractField, field = 1, delim = "-")))
    }
    rownames(counts) <- make.unique(names = as.character(sapply(gene.names, ExtractField, field = 2, delim = "\\t")))
    colnames(counts) <- cell.names

  } else {
    counts <- as.matrix(read.table(matrix.dir, sep = "\t", header = T, row.names = 1))
  }
  colnames(counts) <- make.names(colnames(counts))

  if (make.sparse) {
    counts <- as(counts, "dgCMatrix")
  }
  return(counts)
}


#' Internal function for winsorizing a matrix
#' Adapted from pagoda2: https://github.com/hms-dbmi/pagoda2
winsorize_matrix <- function(mat, trim) {
  if(trim < 1) {
    trim <- trim*ncol(mat)
  }

  mat <- t(mat)
  inplaceWinsorizeSparseCols(mat, trim)
  return(t(mat))
}


#' Filter data by minimum cells, minimum genes, and winsorizes data
#'
#' @param x Input matrix
#' @param min.samples.frac Minimum number of samples for a feature to be included
#' @param trim Fraction of samples to winsorize if trim < 1. If trim > 1, absolute number of samples to winsorize
#' @param min.nonzero.features Minimum number of nonzero features for a sample to be included
#' @param max.sample.sum Samples with greater than this sum of features will be excluded
#'
#' @return Filtered and winsorized data matrix
#'
#' @import Matrix
#' @export
#'
FilterData <- function(x, min.samples.frac, trim, min.nonzero.features = 500, max.sample.sum = 50000) {
  if (is.data.frame(x) | is.matrix(x)) {
    x <- as(as.matrix(x), "dgCMatrix")
  }

  min.cells <- round(ncol(x)*min.samples.frac)
  x <- x[ , Matrix::colSums(x) < max.sample.sum]
  x <- x[ , Matrix::colSums(x > 0) > min.nonzero.features]
  x <- x[Matrix::rowSums(x > 0) > min.cells, ]
  if (trim > 0) {
    x <- t(winsorize_matrix(t(x), trim = trim))
  }
  return(x)
}



#' Normalization and batch effect removal adapted from pagoda2 [insert citation here]
#'
#' @param counts Input data matrix
#' @param depthScale Rescale all columns to this value
#' @param batch Factor the same length as ncol(counts) which contains batch effects to remove
#'
#' @return Normalized data matrix with batch effects removed
#'
#' @import Matrix
#' @export
#'
NormalizeCounts <- function(counts, depthScale = 1e3, batch = NULL) {
  if (is.data.frame(counts) | is.matrix(counts)) {
    counts <- as(as.matrix(counts), "dgCMatrix")
  }

  if(!is.null(batch)) {
    if(!all(colnames(counts) %in% names(batch))) {
      stop("the supplied batch vector doesn't contain all the cells in its names attribute")
    }
    batch <- as.factor(batch[colnames(counts)])
  }

  depth <- Matrix::colSums(counts)
  counts <- Matrix::t(counts)

  if(!is.null(batch)) {
    # dataset-wide gene average
    gene.av <- (Matrix::colSums(counts)+length(levels(batch)))/(sum(depth)+length(levels(batch)))

    # pooled counts, df for all genes
    tc <- colSumByFac(counts, as.integer(batch))[-1,,drop = F]
    tc <- t(log(tc + 1) - log(as.numeric(tapply(depth, batch, sum)) + 1))
    bc <- exp(tc - log(gene.av))

    # adjust every non-0 entry
    count.gene <- rep(1:counts@Dim[2], diff(counts@p))
    counts@x <- counts@x/bc[cbind(count.gene,as.integer(batch)[counts@i + 1])]
  }

  counts <- counts/(depth/depthScale)
  return(Matrix::t(counts))
}


#' Adjust feature variance to remove effect of varying feature means.
#' Adapted from pagoda2: https://github.com/hms-dbmi/pagoda2
#'
#' @param counts Input data matrix
#' @param gam.k Number of additive models to use when fitting mean variance relationship
#' @param plot Whether or not to plot the mean variance relationship
#' @param max.adjusted.variance Maximum adjusted feature variance
#' @param min.adjusted.variance Minimum adjusted feature variance
#' @param verbose Whether or not to print output
#' @param q.val Maximum adjusted p-value for a feature to qualify as overdispersed
#'
#' @return Dataframe with adjusted feature variances
#'
#' @importFrom mgcv gam
#' @export
#'
AdjustVariance <- function(counts, gam.k = 10, plot = F, max.adjusted.variance = 1e3, min.adjusted.variance = 1e-3,
                           verbose = T, q.val = 0.05) {
  if (is.data.frame(counts) | is.matrix(counts)) {
    counts <- as(as.matrix(counts), "dgCMatrix")
  }
  counts <- Matrix::t(counts)

  if(verbose) cat("calculating variance fit ...")
  df <- colMeanVarS(counts, NULL)

  df$m <- log(df$m); df$v <- log(df$v);
  rownames(df) <- colnames(counts);
  vi <- which(is.finite(df$v) & df$nobs >= 0);
  if(length(vi) < gam.k*1.5) { gam.k = 1 };# too few genes
  if(gam.k < 2) {
    if(verbose) cat(" using lm ")
    m <- lm(v ~ m, data = df[vi,])
  } else {
    if(verbose) cat(" using gam ")
    m <- mgcv::gam(as.formula(paste0("v ~ s(m, k = ", gam.k, ")")), data = df[vi,])
  }
  df$res <- -Inf;  df$res[vi] <- resid(m,type='response')
  n.obs <- df$nobs; #diff(counts@p)
  df$lp <- as.numeric(pf(exp(df$res),n.obs,n.obs,lower.tail=F,log.p=F))
  df$lpa <- p.adjust(df$lp, method = "BH")
  n.cells <- nrow(counts)
  df$qv <- as.numeric(qchisq(df$lp, n.cells-1, lower.tail = FALSE,log.p=F)/n.cells)

  ods <- which(df$lpa < q.val)

  df$gsf <- geneScaleFactors <- sqrt(pmax(min.adjusted.variance,pmin(max.adjusted.variance,df$qv))/exp(df$v));
  df$gsf[!is.finite(df$gsf)] <- 0;
  odgenes <- rownames(df)[ods]
  df$overdispersed <- rownames(df) %in% odgenes

  if(plot) {
    smoothScatter(df$m, df$v, main = '', xlab = 'log10[ magnitude ]', ylab = 'log10[ variance ]')
    grid <- seq(min(df$m[vi]), max(df$m[vi]), length.out = 1000)
    lines(grid, predict(m, newdata = data.frame(m = grid)), col = "blue")
    if(length(ods) > 0) {
      points(df$m[ods], df$v[ods], pch = '.', col = 2, cex = 1)
    }
    smoothScatter(df$m[vi], df$qv[vi], xlab = 'log10[ magnitude ]', ylab = '', main = 'adjusted')
    abline(h = 1, lty = 2, col = 8)
    if(is.finite(max.adjusted.variance)) { abline(h = max.adjusted.variance, lty = 2, col = 1) }
    points(df$m[ods], df$qv[ods], col = 2, pch = '.')
  }

  return(df)
}


#' Freeman-Tukey transform for variance stabilization
ft_transform <- function(A) {
  return(sqrt(A) + sqrt(A + 1))
}


#' Normalize, adjust feature variance, and scale data matrix
#'
#' @param counts Input data matrix
#' @param batch Factor with batch identity for each sample
#' @param method Scaling method
#' @param adj.var Whether or not to apply mean variance adjustment for features
#' @param plot.var.adj Whether or not to plot the mean variance relationship
#' @param gam.k Number of additive models to use for variance modeling
#'
#' @return Normalized and scaled data matrix
#'
#' @export
#'
ScaleCounts <- function(counts, batch = NULL, method = "log", adj.var = T, plot.var.adj = F, gam.k = 10) {
  stopifnot(method %in% c("log", "logrank", "ft", "none"))

  scale.factor <- median(Matrix::colSums(counts))
  norm.counts <- NormalizeCounts(counts, batch = batch, depthScale = scale.factor)

  if (adj.var) {
    varinfo.df <- AdjustVariance(norm.counts, plot = plot.var.adj, gam.k = gam.k)
    norm.counts <- norm.counts[rownames(varinfo.df),] * varinfo.df$gsf
  }

  if (method == "logrank") {
    rank.scale <- median(Matrix::colSums(counts > 0))
    norm.counts <- apply(norm.counts, 2, function(x) {
      nz <- which(x > 0)
      x[nz] <- rank(x[nz], ties.method = "average")
      x/length(x) * rank.scale
    })
  } else if (method == "ft") {
    norm.counts@x <- ft_transform(norm.counts@x)
  } else if (method == "log") {
    norm.counts@x <- log(norm.counts@x + 1)
  }

  return(norm.counts)
}


#' Select overdispersed features
#'
#' @param counts Input data matrix
#' @param n.features Number of features to keep
#' @param gam.k Number of additive models to use for variance modeling
#'
#' @import Matrix
#' @export
#'
SelectFeatures <- function(counts, batch = NULL, n.features = 3e3, gam.k = 10) {
  scale.factor <- median(Matrix::colSums(counts))
  norm.counts <- NormalizeCounts(counts, batch = batch, depthScale = scale.factor)

  varinfo <- AdjustVariance(norm.counts, plot = F, gam.k = gam.k)
  varinfo <- varinfo[order(varinfo$lp),]

  return(rownames(varinfo[1:n.features,]))
}


#' Reconstructs the gene expression matrix from the CCA gene loadings and cell embeddings
#'
#' @param se.obj Seurat object with multiple batches aligned via CCA alignment
#' @return Reconstructed gene expression matrix
#'
#' @export
#'
ExtractDebatchedSeurat <- function(se.obj) {
  if (!ncol(se.obj@dr$cca.aligned@cell.embeddings) != length(se.obj@cell.names)) {
    stop("Run CCA alignment first")
  }

  cc.scores <- se.obj@dr$cca.aligned@cell.embeddings
  cc.loadings <- se.obj@dr$cca@gene.loadings[,1:ncol(cc.scores)]

  gene.expr <- cc.loadings %*% t(cc.scores)
  gene.expr <- t(apply(gene.expr, 1, function(x) (x - min(x))/(max(x) - min(x))))

  return(gene.expr)
}


#' Extracts scaled counts from Seurat or Pagoda2 objects
#'
#' @param obj Seurat or Pagoda2 object
#' @param obj.type Must be either "seurat" or "pagoda2"
#' @param rescale Rescale data with plain pagoda2 model
#' @param rescale.method Rescale method, either "log", or "ft"
#' @param batch If rescaling and multiple batches, specify batch
#' @return Normalized gene expression matrix
#'
#' @export
#'
ExtractNormCounts <- function(obj, obj.type = "seurat", rescale = T, rescale.method = "log",
                              batch = NULL) {
  if (!obj.type %in% c("seurat", "pagoda2")) {
    stop("obj.type must be either 'seurat' or 'pagoda2'")
  }

  if (!rescale.method %in% c("log", "ft")) {
    stop("rescale.method must be either 'log', or 'ft'")
  }

  if (obj.type == "seurat") {
    if (rescale) {
      counts <- obj@raw.data[,obj@cell.names]
      norm.counts <- ScaleCounts(counts, batch = batch, method = rescale.method)
    } else {
      norm.counts <- obj@scale.data
      norm.counts <- t(apply(norm.counts, 1, function(x) (x - min(x))/(max(x) - min(x))))
    }
  } else if (obj.type == "pagoda2") {
    if (!is.null(obj$counts)) {
      counts <- Matrix::t(obj$misc$rawCounts[rownames(obj$counts),])
    } else {
      counts <- Matrix::t(obj$misc$rawCounts)
    }
    norm.counts <- ScaleCounts(counts, batch = batch, method = rescale.method)
  }

  return(norm.counts)
}
