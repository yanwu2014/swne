## Matrix processing functions

.extract_field <- function (string, field = 1, delim = "_") {
  fields <- as.numeric(unlist(strsplit(x = as.character(x = field), split = ",")))
  if (length(fields) == 1) {
    return(strsplit(string, split = delim)[[1]][field])
  }
  return(paste(strsplit(string, split = delim)[[1]][fields], collapse = delim))
}


ReadData <- function(matrix.dir, make.sparse = T) {
  if (dir.exists(matrix.dir)) {
    if (!grepl("\\/$", matrix.dir)) { matrix.dir <- paste(matrix.dir, "/", sep = "") }

    barcode.loc <- paste0(matrix.dir, "barcodes.tsv")
    gene.loc <- paste0(matrix.dir, "genes.tsv")
    matrix.loc <- paste0(matrix.dir, "matrix.mtx")

    if (!file.exists(barcode.loc)) { stop("Barcode file missing") }
    if (!file.exists(gene.loc)) { stop("Gene name file missing") }
    if (!file.exists(matrix.loc)) { stop("Expression matrix file missing") }

    counts <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    gene.names <- readLines(gene.loc)
    if (all(grepl(pattern = "\\-1$", cell.names))) {
      cell.names <- as.vector(as.character(sapply(cell.names, .extract_field, field = 1, delim = "-")))
    }
    rownames(counts) <- make.unique(names = as.character(sapply(gene.names, .extract_field, field = 2, delim = "\\t")))
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


.winsorize_matrix <- function(mat, trim) {
  if(trim < 1) {
    trim <- trim*ncol(mat)
  }

  mat <- t(mat)
  inplaceWinsorizeSparseCols(mat, trim)
  return(t(mat))
}


# Filter data by minimum cells, minimum genes
FilterData <- function(x, min.cells.frac, trim, min.genes = 500, min.expr = 0, max.transcripts = 70000) {
  min.cells <- round(ncol(x)*min.cells.frac)
  x <- x[ , Matrix::colSums(x) < max.transcripts]
  x <- x[ , Matrix::colSums(x > min.expr) > min.genes]
  x <- x[Matrix::rowSums(x > min.expr) > min.cells, ]
  if (trim > 0) {
    x <- t(.winsorize_matrix(t(x), trim = trim))
  }
  return(x)
}



# Normalization and batch effect removal
NormalizeCounts <- function(counts, depthScale = 1e3, batch = NULL) {
  if(!is.null(batch)) {
    if(!all(colnames(counts) %in% names(batch))) {
      stop("the supplied batch vector doesn't contain all the cells in its names attribute")
    }
    batch <- as.factor(batch[colnames(counts)])
  }

  depth <- Matrix::colSums(counts)
  counts <- t(counts)

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
  return(t(counts))
}



AdjustVariance <- function(counts, gam.k = 5, alpha = 5e-2, plot = F, use.unadjusted.pvals = F, max.adjusted.variance = 1e3,
                           min.adjusted.variance = 1e-3, verbose = T) {
  counts <- t(counts)

  if(verbose) cat("calculating variance fit ...")
  df <- colMeanVarS(counts, NULL)

  df$m <- log(df$m); df$v <- log(df$v);
  rownames(df) <- colnames(counts);
  vi <- which(is.finite(df$v) & df$nobs >= 0);
  if(length(vi) < gam.k*1.5) { gam.k=1 };# too few genes
  if(gam.k < 2) {
    if(verbose) cat(" using lm ")
    m <- lm(v ~ m, data = df[vi,])
  } else {
    if(verbose) cat(" using gam ")
    require(mgcv)
    m <- mgcv::gam(v ~ s(m, k = gam.k), data = df[vi,])
  }
  df$res <- -Inf;  df$res[vi] <- resid(m,type='response')
  n.obs <- df$nobs; #diff(counts@p)
  df$lp <- as.numeric(pf(exp(df$res),n.obs,n.obs,lower.tail=F,log.p=F))
  df$lpa <- p.adjust(df$lp, method = "BH")
  n.cells <- nrow(counts)
  df$qv <- as.numeric(qchisq(df$lp, n.cells-1, lower.tail = FALSE,log.p=F)/n.cells)

  if(use.unadjusted.pvals) {
    ods <- which(df$lp < alpha)
  } else {
    ods <- which(df$lpa < alpha)
  }

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


.ft_transform <- function(A) {
  return(sqrt(A) + sqrt(A + 1))
}


ScaleCounts <- function(counts, batch = NULL, method = "log", adj.var = T, plot.var.adj = F) {
  stopifnot(method %in% c("log", "logrank", "ft", "none"))

  scale.factor <- median(Matrix::colSums(counts))
  norm.counts <- NormalizeCounts(counts, batch = batch, depthScale = scale.factor)

  if (adj.var) {
    varinfo.df <- AdjustVariance(norm.counts, plot = plot.var.adj)
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
