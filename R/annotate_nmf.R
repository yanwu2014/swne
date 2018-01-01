## Geneset handling functions

ReadGenesets <- function(geneset.name) {
  geneset.to.use <- as.list(readLines(geneset.name))
  geneset.to.use <- lapply(geneset.to.use, function (v) strsplit(v, '\t')[[1]])
  geneset.names <- unlist(lapply(geneset.to.use, function(x) x[[1]]))
  geneset.to.use <- lapply(geneset.to.use, function(v) v[3:length(v)])
  names(geneset.to.use) <- geneset.names
  names(geneset.to.use) <- sapply(names(geneset.to.use), function(x) gsub(" ", "_", x))
  return(geneset.to.use)
}


## Helper function for filtering genesets
.clean_genesets <- function(go.env, min.size = 5, max.size = 500, annot = FALSE) {
  go.env <- as.list(go.env)
  size <- unlist(lapply(go.env, length))
  go.env <- go.env[size > min.size & size < max.size]
  return(go.env)
}


FilterGenesets <- function(geneset.to.use, gene.names, min.size = 5, max.size = 500) {
  geneset.to.use <- lapply(geneset.to.use, function (x) return(x[x %in% gene.names]))
  geneset.to.use <- .clean_genesets(geneset.to.use, min.size = min.size, max.size = max.size)
  return(geneset.to.use)
}


WriteGenesets <- function(genesets, file.name) {
  genesets <- lapply(names(genesets), function(name) {x <- genesets[[name]]; x <- c(name,name,x); return(x);})
  n.cols <- 1.5*max(unlist(lapply(genesets, length)))
  empty <- lapply(genesets, write, file.name, append=TRUE, ncolumns = n.cols, sep = '\t')
}


## Functions for annotating NMF components

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


ProjectGenesets <- function(norm.counts, genesets, loss = "mse") {
  full.genesets.matrix <- .genesets_indicator(genesets, inv = F, return.numeric = T)
  project_nmf(norm.counts[rownames(full.genesets.matrix),], full.genesets.matrix, loss = loss)
}


ComponentAssociation <- function(X, Y, n.cores = 8, metric = "IC") {
  cl <- snow::makeCluster(n.cores, type = "SOCK")
  snow::clusterExport(cl, c("Y", "MutualInf"), envir = environment())
  if (metric == "IC") {
    source.log <- snow::parLapply(cl, 1:length(cl), function(i) library(MASS))
    assoc <- t(snow::parApply(cl, X, 1, function(v) apply(Y, 1, function(u) MutualInf(u, v))))
  } else if (metric == "pearson") {
    assoc <- t(snow::parApply(cl, X, 1, function(v) apply(Y, 1, function(u) cor(u, v))))
  } else if (metric == "spearman") {
    assoc <- t(snow::parApply(cl, X, 1, function(v) apply(Y, 1, function(u) cor(u, v, method = "spearman"))))
  } else {
    stop("Invalid correlation metric")
  }
  stopCluster(cl)

  rownames(assoc) <- rownames(X)
  colnames(assoc) <- rownames(Y)
  return(assoc)
}



SummarizeAssocGenes <- function(gene.nmf.assoc, genes.return = 10, genes.use = NULL) {
  if (!is.null(genes.use)) {
    gene.nmf.assoc <- gene.nmf.assoc[genes.use,]
  }

  nmf.genes.df <- do.call("rbind", lapply(1:ncol(gene.nmf.assoc), function(i) {
    genes.df <- data.frame(Cor = gene.nmf.assoc[,i])
    genes.df$gene <- rownames(gene.nmf.assoc)
    genes.df$Metagene <- colnames(gene.nmf.assoc)[[i]]
    genes.df <- genes.df[order(genes.df$Cor, decreasing = T),]
    head(genes.df, n = genes.return)
  }))

  rownames(nmf.genes.df) <- NULL
  return(nmf.genes.df)
}


SummarizeAssocGenesets <- function(genesets.nmf.assoc, n.return = 5) {
  assoc.genesets.df <- do.call("rbind", lapply(colnames(genesets.nmf.assoc), function(nf) {
    df <- data.frame(Cor = genesets.nmf.assoc[ ,nf])
    df$Genesets <- rownames(df)
    df$NMF <- nf
    head(df[order(df$Cor, decreasing = T), ], n = n.return)
  })); rownames(assoc.genesets.df) <- NULL;
  assoc.genesets.df
}


# Compute Information Coefficient [IC]
# Pablo Tamayo Dec 30, 2015
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
    kde2d.xy <- kde2d(x, y, n = n.grid, h = delta)
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

FlattenGroups <- function(groups.list) {
  cell.names <- unlist(groups.list, F, F)
  if (length(cell.names) > length(unique(cell.names))) {
    print("Warning: removing all cells with duals")
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


UnflattenGroups <- function(groups, min.cells = 1) {
  groups.list <- c()
  unique.groups <- unique(groups)
  for (g in unique.groups) {
    g.cells <- names(groups[groups == g])
    if (length(g.cells) > min.cells) {
      groups.list[[g]] <- g.cells
    }
  }
  return(groups.list)
}


ReadGroups <- function(groups.file, sep.char = ",") {
  group.data <- read.table(groups.file, sep = sep.char, header = F, stringsAsFactors = F)
  groups <- group.data[[1]]
  full.groups.list <- lapply(rownames(group.data), function(i) sapply(strsplit(group.data[i,2], split = ',')[[1]], trimws))
  full.groups.list <- lapply(full.groups.list, function(x) make.names(x))
  names(full.groups.list) <- groups
  return(full.groups.list)
}


WriteGroups <- function(groups.list, out.file) {
  group.data <- sapply(groups.list, function(x) paste('\"', paste(x, collapse = ", "), '\"', sep = ""))
  group.data <- sapply(names(group.data), function(x) paste(x, group.data[[x]], sep = ","))
  fileConn = file(out.file)
  writeLines(group.data, fileConn)
  close(fileConn)
}
