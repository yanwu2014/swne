#' @import methods
#' @import ggplot2
#' @import RColorBrewer
#' @import Matrix
#' @import NNLM
#' @import ggrepel
#' @import ggConvexHull
#' @import snow
#' @useDynLib swne
#' @importFrom Rcpp evalCpp
NULL

#### SWNE functions ####

.normalize_vector <- function(x, method = "scale", n_ranks = 10000) {
  stopifnot(method %in% c("rank", "scale", "bounded"))
  if (method == "scale") {
    x.scale <- (x - mean(x))/sd(x)
  } else if (method == "rank") {
    x.scale <- rank(x) / length(x) * n_ranks
  } else if (method == "bounded") {
    x.max <- max(x)
    x.min <- min(x)
    x.scale <- (x - x.min)/(x.max - x.min)
  }

  return(x.scale)
}



.get_component_coords <- function(H, distance = "pearson") {
  H <- t(H)
  stopifnot(distance %in% c("pearson", "IC"))
  if (distance == "pearson") {
    H.dist <- sqrt(2*(1 - cor(H)))
  } else if (distance == "IC") {
    require(usedist)
    H.dist <- sqrt(2*(1 - as.matrix(dist_make(t(H), distance_fcn = MutualInf, method = "IC"))))
  }

  H.coords <- MASS::sammon(H.dist, k = 2, niter = 250)$points
  H.coords <- apply(H.coords, 2, .normalize_vector, method = "bounded")
  colnames(H.coords) <- c("x","y")

  return(H.coords)
}



.get_sample_coords <- function(H, H.coords, alpha = 1, n_pull = NULL) {
  if (n_pull < 3) { n_pull <- 3 }
  if (is.null(n_pull) || n_pull > nrow(H.coords)) { n_pull <- nrow(H.coords) }

  sample.coords <- t(sapply(1:ncol(H), function(i) {
    pull.ix <- order(H[,i], decreasing = T)[1:n_pull]
    pull.sum <- sum(H[pull.ix,i]^alpha)
    x <- sum(sapply(pull.ix, function(j) {
      H.coords[j,1] * ( H[j,i]^alpha/pull.sum )
    }))
    y <- sum(sapply(pull.ix, function(j) {
      H.coords[j,2] * ( H[j,i]^alpha/pull.sum )
    }))
    return(c(x,y))
  }))
  colnames(sample.coords) <- c("x","y")
  rownames(sample.coords) <- colnames(H)
  return(sample.coords)
}


#' Embeds NMF components and samples in a 2D space
#'
#' @param H NMF scores (NMF components x samples)
#' @param SNN Shared nearest neighbors matrix (or other similarity matrix)
#' @param alpha.exp Increasing alpha.exp increases how much the NMF components "pull" the samples
#' @param snn.exp Decreasing snn.exp increases the effect of the similarity matrix on the embedding
#' @param n_pull Number of NMF components pulling on each sample. n_pull >= 3
#' @param dist.use Similarity function to use for calculating NMF component positions. Options include pearson and IC.
#' @param min.snn Minimum SNN value
#'
#' @return A list of nmf component coordinates in 2D space (H.coords) and sample coordinates (sample.coords)
#'
#' @export
#'
EmbedSWNE <- function(H, SNN = NULL, alpha.exp = 1, snn.exp = 1.0, n_pull = NULL, dist.use = "pearson",
                      min.snn = 0.05) {
  H <- H[ ,colSums(H) > 0]
  H.coords <- .get_component_coords(H, distance = dist.use)
  H.coords <- data.frame(H.coords); H.coords$name <- rownames(H.coords);

  sample.coords <- .get_sample_coords(H, H.coords, alpha = alpha.exp, n_pull = n_pull)

  if (!is.null(SNN)) {
    SNN@x[SNN@x < min.snn] <- 0
    SNN@x <- SNN@x^snn.exp
    SNN <- SNN/Matrix::rowSums(SNN)

    x <- sapply(1:nrow(SNN), function(i) sum(SNN[i,]*sample.coords[,1]))
    y <- sapply(1:nrow(SNN), function(i) sum(SNN[i,]*sample.coords[,2]))
    sample.coords <- as.matrix(data.frame(x, y)); rownames(sample.coords) <- colnames(H);
  }

  return(list(H.coords = H.coords, sample.coords = data.frame(sample.coords)))
}


#' Projects new data onto an existing swne embedding
#'
#' @param swne.embedding Existing swne embedding from EmbedSWNE
#' @param H.test Projected test data NMF scores
#' @param SNN.test Projected test data shared nearest neighbors matrix
#' @param alpha.exp Increasing alpha.exp increases how much the NMF components "pull" the samples
#' @param snn.exp Decreasing snn.exp increases the effect of the similarity matrix on the embedding
#' @param n_pull Number of NMF components pulling on each sample. n_pull >= 3
#'
#' @return Matrix of sample embeddings
#'
#' @export
#'
ProjectSWNE <- function(swne.embedding, H.test, SNN.test = NULL, alpha.exp = 1, snn.exp = 1, n_pull = NULL) {
  sample.coords.test <- .get_sample_coords(H.test, swne.embedding$H.coords, alpha = alpha.exp, n_pull = n_pull)
  sample.coords.train <- as.matrix(swne.embedding$sample.coords)
  sample.coords.all <- rbind(sample.coords.test, sample.coords.train)

  if (!is.null(SNN.test)) {
    SNN.test <- SNN.test[rownames(sample.coords.test), rownames(sample.coords.train)]
    SNN.diag <- Matrix::Diagonal(ncol(H.test)); rownames(SNN.diag) <- colnames(SNN.diag) <- colnames(H.test);
    SNN.test <- cbind(SNN.diag, SNN.test); rm(SNN.diag);
    SNN.test <- SNN.test^snn.exp
    SNN.test <- SNN.test/Matrix::rowSums(SNN.test)

    x <- sapply(1:nrow(SNN.test), function(i) sum(SNN.test[i,]*sample.coords.all[,1]))
    y <- sapply(1:nrow(SNN.test), function(i) sum(SNN.test[i,]*sample.coords.all[,2]))
    sample.coords.test <- data.frame(x, y); rownames(sample.coords.test) <- colnames(H.test);
  } else {
    sample.coords.test <- data.frame(sample.coords.test)
  }
  return(sample.coords.test)
}


#' Plots swne embedding
#'
#' @param swne.embedding SWNE embedding (list of NMF and sample coordinates) from EmbedSWNE
#' @param alpha.plot Data point transparency
#' @param sample.groups Factor defining sample groups
#' @param do.label Label the sample groups
#' @param label.size Label font size
#' @param pt.size Sample point size
#' @param samples.plot Vector of samples to plot. Default is NULL, which plots all samples.
#' @param show.legend If sample groups defined, show legend
#' @param seed Seed for sample group color reproducibility
#'
#' @return ggplot2 object with swne plot
#'
#' @export
#'
PlotSWNE <- function(swne.embedding, alpha.plot = 0.25, sample.groups = NULL, do.label = F,
                     label.size = 4.5, pt.size = 1, samples.plot = NULL, show.legend = T,
                     seed = NULL) {
  H.coords <- swne.embedding$H.coords
  sample.coords <- swne.embedding$sample.coords

  sample.groups <- factor(sample.groups[rownames(sample.coords)])
  sample.coords$pt.size <- pt.size

  if (!is.null(sample.groups)) {
    set.seed(seed)
    sample.groups <- factor(sample.groups, levels = sample(levels(sample.groups)))
    sample.coords$sample.groups <- sample.groups[rownames(sample.coords)]
  } else {
    sample.coords$sample.groups <- factor(rep(1, nrow(sample.coords)))
  }

  if (!is.null(samples.plot)) {
    sample.coords <- sample.coords[samples.plot,]
    sample.coords <- droplevels(sample.coords)
  }

  ## Plot sample coordinates
  ggobj <- ggplot() +
    geom_point(data = sample.coords, aes(x, y, colour = sample.groups, fill = sample.groups),
               alpha = alpha.plot, size = pt.size) +
    theme_classic() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                            axis.ticks = element_blank(), axis.line = element_blank(),
                            axis.text = element_blank()) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = label.size)))

  ## Plot NMF points and draw convex hull
  ggobj <- ggobj +
    ggConvexHull::geom_convexhull(data = H.coords, aes(x, y), alpha = 0.1, fill = NA, size = 0.75,
                                  colour = "grey", linetype = "dotdash") +
    geom_point(data = H.coords, aes(x, y), size = 3, color = "blue")

  ## Plot text labels
  if (do.label && is.factor(sample.groups)) {
    group.pts.x <- tapply(sample.coords$x, sample.coords$sample.groups, median)
    group.pts.y <- tapply(sample.coords$y, sample.coords$sample.groups, median)

    group.pts <- data.frame(x = group.pts.x, y = group.pts.y, name = levels(sample.coords$sample.groups))
    label.pts <- rbind(H.coords, group.pts)

    ggobj <- ggobj + ggrepel::geom_text_repel(data = label.pts, mapping = aes(x, y, label = name),
                                              size = label.size, box.padding = 0.15)
  } else {
    ggobj <- ggobj + ggrepel::geom_text_repel(data = H.coords, mapping = aes(x, y, label = name),
                                              size = label.size, box.padding = 0.15)
  }

  if (!show.legend) {
    ggobj <- ggobj + theme(legend.position = "none")
  }

  return(ggobj)
}


#' Plots swne embedding with feature overlayed
#'
#' @param H.coords NMF coordinates
#' @param sample.coords Sample coordinates
#' @param feature.score Feature vector to overlay
#' @param n.colors Number of colors to bin feature vector into
#' @param alpha.plot Data point transparency
#' @param quantiles Quantiles to trim outliers from
#' @param samples.plot Samples to actually plot. Default is NULL, which plots all samples
#' @param label.size Label font size
#' @param pt.size Sample point size
#'
#' @return ggplot2 object with swne plot with feature overlayed
#'
#' @export
#'
FeaturePlotSWNE <- function(H.coords, sample.coords, feature.scores, n.colors = 5, alpha.plot = 0.5,
                            quantiles = c(0.05, 0.95), samples.plot = NULL, label.size = 4.5, pt.size = 1) {
  feature.scores <- as.numeric(feature.scores[rownames(sample.coords)])
  feature.quantiles <- quantile(feature.scores, probs = quantiles)

  min.quantile <- feature.quantiles[[1]]; max.quantile <- feature.quantiles[[2]];
  feature.scores[feature.scores < min.quantile] <- min.quantile
  feature.scores[feature.scores > max.quantile] <- max.quantile
  feature.cut <- as.numeric(as.factor(cut(feature.scores, breaks = n.colors)))

  sample.coords$pt.size <- pt.size
  sample.coords$feature <- feature.scores
  sample.coords$sample.groups <- as.factor(feature.cut)

  if (!is.null(samples.plot)) {
    sample.coords <- sample.coords[samples.plot,]
    sample.coords <- droplevels(sample.coords)
  }

  ## Plot sample coordinates
  ggobj <- ggplot() +
    geom_point(data = sample.coords, aes(x, y, colour = sample.groups, fill = sample.groups), alpha = alpha.plot, size = pt.size) +
    theme_classic() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                            axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank()) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = label.size))) +
    scale_colour_brewer(palette = "Blues")

  ## Plot NMF points and draw convex hull
  ggobj <- ggobj +
    ggConvexHull::geom_convexhull(data = H.coords, aes(x, y), alpha = 0.1, fill = NA, size = 0.75,
                                  colour = "grey", linetype = "dotdash") +
    geom_point(data = H.coords, aes(x, y), size = 3, color = "blue")

  ## Plot text labels
  ggobj <- ggobj + ggrepel::geom_text_repel(data = H.coords, mapping = aes(x, y, label = name),
                                            size = label.size, box.padding = 0.15)
  ggobj <- ggobj + theme(legend.position = "none")

  return(ggobj)
}


#' Plots 2D embedding
#'
#' @param dim.scores 2D embedding coordinates. Must be 2 x n.samples
#' @param sample.groups Factor defining sample groups
#' @param x.lab X axis label
#' @param y.lab Y axis label
#' @param main.title Main plot title
#' @param pt.size Sample point size
#' @param font.size Font size for axis labels
#' @param alpha.plot Data point transparency
#' @param do.label Label the sample groups
#' @param label.size Label font size
#' @param show.legend If sample groups defined, show legend
#' @param show.axes Plot x and y axes
#' @param seed Seed for sample group color reproducibility
#'
#' @return ggplot2 object with 2d plot
#'
#' @export
#'
PlotDims <- function(dim.scores, sample.groups = NULL, x.lab = "tsne1", y.lab = "tsne2",
                     main.title = NULL, pt.size = 1.0, font.size = 12, alpha.plot = 1.0, do.label = T,
                     label.size = 4, show.legend = T, show.axes = T, seed = NULL) {
  if (!is.null(sample.groups)) {
    set.seed(seed)
    sample.groups <- factor(sample.groups, levels = sample(levels(sample.groups)))
  }

  gg.df <- data.frame(x = dim.scores[,1], y = dim.scores[,2], sample.groups = sample.groups)
  ggobj <- ggplot(gg.df, aes(x, y)) + geom_point(size = pt.size, alpha = alpha.plot, aes(colour = sample.groups)) +
    theme(text = element_text(size = font.size)) +
    xlab(x.lab) + ylab(y.lab) + ggtitle(main.title) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = label.size)))

  if (!show.axes) { ggobj <- ggobj + theme_void() }

  if (do.label && is.factor(sample.groups)) {
    group.pts.x <- tapply(gg.df$x, gg.df$sample.groups, function(v) median(v))
    group.pts.y <- tapply(gg.df$y, gg.df$sample.groups, function(v) median(v))
    group.pts <- data.frame(x = group.pts.x, y = group.pts.y)
    group.pts$ident <- levels(gg.df$sample.groups)
  }

  ggobj <- ggobj + ggrepel::geom_text_repel(data = group.pts, mapping = aes(label = ident), size = label.size,
                                            box.padding = 0.15)

  if (!show.legend) { ggobj <- ggobj + theme(legend.position = "none") }

  ggobj
}


#' Plots a heatmap using ggplot2. Adapted from [site]
#'
#' @param m Matrix to plot heatmap for
#' @param rescaling Scale by row or columns
#' @param clustering Cluster by row or columns
#' @param labCol Label columns
#' @param labRow Label rows
#' @param border Plot heatmap border
#' @param heatscale Color map
#' @param legend.title Title of heatmap legend
#' @param x.lab.size x label size
#' @param y.lab.size y label size
#'
#' @return ggplot2 heatmap
#'
#' @import reshape
#'
#' @export
#'
ggHeat <- function(m, rescaling = 'none', clustering = 'none', labCol = T, labRow = T, border = F,
                   heatscale = c(low = 'skyblue', mid = 'white', high = 'tomato'), legend.title = NULL, x.lab.size = 10,
                   y.lab.size = 10) {
  ## you can either scale by row or column not both!
  ## if you wish to scale by both or use a differen scale method then simply supply a scale
  ## function instead NB scale is a base funct

  if(is.function(rescaling)) {
    m = rescaling(m)
  } else
  {
    if(rescaling == 'column')
      m = scale(m, center = T)
    if(rescaling == 'row')
      m = t(scale(t(m), center = T))
  }

  ## I have supplied the default cluster and euclidean distance- and chose to cluster after scaling
  ## if you want a different distance/cluster method-- or to cluster and then scale
  ## then you can supply a custom function

  if(is.function(clustering)) {
    m = clustering(m)
  } else {
    if(clustering == 'row')
      m = m[hclust(dist(m))$order, ]
    if(clustering == 'column')
      m = m[ ,hclust(dist(t(m)))$order]
    if(clustering=='both')
      m = m[hclust(dist(m))$order, hclust(dist(t(m)))$order]
  }
  ## this is just reshaping into a ggplot format matrix and making a ggplot layer

  rows = dim(m)[1]
  cols = dim(m)[2]
  melt.m = cbind(rowInd = rep(1:rows, times = cols), colInd = rep(1:cols, each = rows), reshape::melt(m))
  g = ggplot(data = melt.m)

  ## add the heat tiles with or without a white border for clarity
  if(border == TRUE)
    g2 = g + geom_rect(aes(xmin = colInd - 1, xmax = colInd, ymin = rowInd - 1, ymax = rowInd, fill = value), colour = 'white')
  if(border == FALSE)
    g2 = g + geom_rect(aes(xmin = colInd - 1, xmax = colInd, ymin = rowInd - 1, ymax = rowInd, fill = value))

  ## add axis labels either supplied or from the colnames rownames of the matrix
  if(labCol == T)
    g2 = g2 + scale_x_continuous(breaks = (1:cols) - 0.5, labels = colnames(m), expand = c(0.005,0))
  if(labCol == F)
    g2 = g2 + scale_x_continuous(breaks = (1:cols) - 0.5, labels = rep('', cols))
  if(labRow == T)
    g2 = g2 + scale_y_continuous(breaks = (1:rows) - 0.5, labels = rownames(m), expand = c(0.005,0))
  if(labRow == F)
    g2 = g2 + scale_y_continuous(breaks = (1:rows) - 0.5, labels = rep('', rows))

  ## get rid of grey panel background and gridlines
  g2 = g2 + theme(panel.grid.minor = element_line(colour = NA), panel.grid.major = element_line(colour=NA),
                  panel.background = element_rect(fill = NA, colour = NA), axis.text.x = element_text(angle = 90, hjust = 1, size = x.lab.size),
                  axis.ticks = element_blank(), axis.text.y = element_text(size = y.lab.size))

  ## finally add the fill colour ramp of your choice (default is blue to red)-- and return
  return(g2 + scale_fill_gradient2(low = heatscale[1], mid = heatscale[2], high = heatscale[3], guide = guide_colorbar(title = legend.title)))
}

