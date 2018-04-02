#' @import methods
#' @import ggplot2
#' @useDynLib swne
#' @importFrom Rcpp evalCpp
NULL

#### SWNE functions ####

#' Helper function for normalizing vectors with different methods.
normalize_vector <- function(x, method = "scale", n_ranks = 10000) {
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


#' Calculates the coordinates of the NMF factors via Sammon mapping
#'
#' @importFrom usedist dist_make
#' @importFrom proxy simil
#'
get_factor_coords <- function(H, distance = "IC") {
  H <- t(H)
  stopifnot(distance %in% c("pearson", "IC", "cosine"))
  if (distance == "pearson") {
    H.dist <- sqrt(2*(1 - cor(H)))
  } else if (distance == "IC") {
    H.dist <- sqrt(2*(1 - as.matrix(usedist::dist_make(t(H), distance_fcn = MutualInf, method = "IC"))))
  } else if (distance == "cosine") {
    H.dist <- sqrt(2*(1 - proxy::simil(H, method = "cosine", by_rows = F)))
  }

  H.coords <- MASS::sammon(H.dist, k = 2, niter = 250)$points
  H.coords <- apply(H.coords, 2, normalize_vector, method = "bounded")
  colnames(H.coords) <- c("x","y")

  return(H.coords)
}


#' Calculate sample coordinates using the NMF scores and NMF factor coordinates
get_sample_coords <- function(H, H.coords, alpha = 1, n_pull = NULL) {
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


#' Projects NMF factors and samples in a 2D
#'
#' @param H NMF factors (factors x samples)
#' @param SNN Shared nearest neighbors matrix (or other similarity matrix)
#' @param alpha.exp Increasing alpha.exp increases how much the NMF factors "pull" the samples
#' @param snn.exp Decreasing snn.exp increases the effect of the similarity matrix on the embedding
#' @param n_pull Number of factors pulling on each sample. Must be >= 3
#' @param dist.use Similarity function to use for calculating factor positions. Options include pearson, IC, cosine.
#' @param min.snn Minimum SNN value
#'
#' @return A list of factor (H.coords) and sample coordinates (sample.coords) in 2D
#'
#' @export
#'
EmbedSWNE <- function(H, SNN = NULL, alpha.exp = 1, snn.exp = 1.0, n_pull = NULL, dist.use = "IC",
                      min.snn = 0.0) {
  H <- H[ ,colSums(H) > 0]
  H.coords <- get_factor_coords(H, distance = dist.use)
  H.coords <- data.frame(H.coords)
  H.coords$name <- rownames(H.coords)
  rownames(H.coords) <- NULL

  sample.coords <- get_sample_coords(H, H.coords, alpha = alpha.exp, n_pull = n_pull)

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


#' Embeds features relative to factor coordinates
#'
#' @param swne.embedding Existing swne embedding from EmbedSWNE
#' @param feature.assoc Feature loadings or correlations (features x factors)
#' @param features.embed Names of features to embed
#' @param alpha.exp Increasing alpha.exp increases how much the factors "pull" the features
#' @param n_pull Number of factors pulling on each feature. Must be >= 3
#' @param scale.cols Whether or not to scale the input columns to 0 - 1
#' @return swne.embedding with feature coordinates (feature.coords)
#'
#' @export
#'
EmbedFeatures <- function(swne.embedding, feature.assoc, features.embed, alpha.exp = 1, n_pull = NULL,
                          scale.cols = T) {
  feature.assoc <- t(feature.assoc[features.embed,])
  stopifnot(nrow(swne.embedding$H.coords) == nrow(feature.assoc))

  if (scale.cols) {
    feature.assoc <- apply(feature.assoc, 2, normalize_vector, method = "bounded")
  }

  feature.coords <- get_sample_coords(feature.assoc, swne.embedding$H.coords, alpha = alpha.exp, n_pull = n_pull)
  feature.coords <- data.frame(feature.coords)
  feature.coords$name <- rownames(feature.coords); rownames(feature.coords) <- NULL;

  swne.embedding$feature.coords <- feature.coords
  return(swne.embedding)
}


#' Projects new data onto an existing swne embedding
#'
#' @param swne.embedding Existing swne embedding from EmbedSWNE
#' @param H.test Test factor scores
#' @param SNN.test Test SNN matrix
#' @param alpha.exp Increasing alpha.exp increases how much the factors "pull" the samples
#' @param snn.exp Decreasing snn.exp increases the effect of the similarity matrix on the embedding
#' @param n_pull Number of factors pulling on each sample. Must be >= 3
#'
#' @return Matrix of sample embeddings
#'
#' @export
#'
ProjectSWNE <- function(swne.embedding, H.test, SNN.test = NULL, alpha.exp = 1, snn.exp = 1, n_pull = NULL) {
  sample.coords.test <- get_sample_coords(H.test, swne.embedding$H.coords, alpha = alpha.exp, n_pull = n_pull)
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


#' Function for renaming NMF factors to something more interpretable.
#' If the NMF factor name is the empty string, "", then the NMF factor will not be plotted
#'
#' @param swne.embedding List of NMF and sample coordinates from EmbedSWNE
#' @param new.names Named character vector with the old values as names and the new values as values
#' @param set.empty If true, any old NMF names that weren't renamed are set to the empty string
#'
#' @return SWNE embedding with NMFs renamed.
#'
#' @importFrom plyr revalue
#' @export
#'
RenameFactors <- function(swne.embedding, name.mapping, set.empty = T) {
  old.names <- swne.embedding$H.coords$name
  new.names <- revalue(old.names, name.mapping)
  if (set.empty) { new.names[new.names == old.names] <- "" }

  swne.embedding$H.coords$name <- new.names
  return(swne.embedding)
}


#' Plots swne embedding
#'
#' @param swne.embedding SWNE embedding (list of factor and sample coordinates) from EmbedSWNE
#' @param alpha.plot Data point transparency
#' @param sample.groups Factor defining sample groups
#' @param do.label Label the sample groups
#' @param label.size Label font size
#' @param pt.size Sample point size
#' @param samples.plot Vector of samples to plot. Default is NULL, which plots all samples.
#' @param show.legend If sample groups defined, show legend
#' @param seed Seed for sample groups color reproducibility
#'
#' @return ggplot2 object with swne plot
#'
#' @import ggrepel
#' @export
#'
PlotSWNE <- function(swne.embedding, alpha.plot = 0.25, sample.groups = NULL, do.label = F,
                     label.size = 4.5, pt.size = 1, samples.plot = NULL, show.legend = T,
                     seed = NULL) {
  H.coords <- swne.embedding$H.coords
  H.coords.plot <- subset(H.coords, name != "")
  sample.coords <- swne.embedding$sample.coords
  feature.coords <- swne.embedding$feature.coords

  sample.groups <- factor(sample.groups[rownames(sample.coords)])
  sample.coords$pt.size <- pt.size

  set.seed(seed)
  if (!is.null(sample.groups)) {
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
                            axis.text = element_blank(), legend.title = element_blank()) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = label.size)))

  ## Plot factors
  if (nrow(H.coords.plot) > 0) {
    ggobj <- ggobj + geom_point(data = H.coords.plot, aes(x, y), size = 2.5, color = "darkblue")
  }

  if (!is.null(feature.coords)) {
    ggobj <- ggobj + geom_point(data = feature.coords, aes(x, y), size = 2.5, color = "darkred")
  }

  ## Plot text labels
  if (do.label && is.factor(sample.groups)) {
    group.pts.x <- tapply(sample.coords$x, sample.coords$sample.groups, median)
    group.pts.y <- tapply(sample.coords$y, sample.coords$sample.groups, median)

    group.pts <- data.frame(x = group.pts.x, y = group.pts.y, name = levels(sample.coords$sample.groups))
    label.pts <- rbind(H.coords.plot, group.pts, feature.coords)

    ggobj <- ggobj + ggrepel::geom_text_repel(data = label.pts, mapping = aes(x, y, label = name),
                                              size = label.size, box.padding = 0.15)
  } else if (nrow(H.coords.plot) > 0) {
    label.pts <- rbind(H.coords.plot, feature.coords)
    ggobj <- ggobj + ggrepel::geom_text_repel(data = label.pts, mapping = aes(x, y, label = name),
                                              size = label.size, box.padding = 0.15)
  } else if (!is.null(feature.coords)) {
    ggobj <- ggobj + ggrepel::geom_text_repel(data = feature.coords, mapping = aes(x, y, label = name),
                                              size = label.size, box.padding = 0.15)
  }

  if (!show.legend) {
    ggobj <- ggobj + theme(legend.position = "none")
  }

  return(ggobj)
}


#' Plots swne embedding with feature overlayed
#'
#' @param swne.embedding SWNE embedding (list of NMF and sample coordinates) from EmbedSWNE
#' @param feature.score Feature vector to overlay
#' @param feature.name Name of feature to overlay
#' @param alpha.plot Data point transparency
#' @param quantiles Quantiles to trim outliers from
#' @param samples.plot Samples to actually plot. Default is NULL, which plots all samples
#' @param label.size Label font size
#' @param pt.size Sample point size
#' @param color.palette RColorbrewer palette to use
#'
#' @return ggplot2 object with swne plot with feature overlayed
#'
#' @export
#'
FeaturePlotSWNE <- function(swne.embedding, feature.scores, feature.name = NULL, alpha.plot = 0.5,
                            quantiles = c(0.01, 0.99), samples.plot = NULL, label.size = 4.5,
                            pt.size = 1, color.palette = "YlOrRd") {
  H.coords <- swne.embedding$H.coords
  H.coords.plot <- subset(H.coords, name != "")
  sample.coords <- swne.embedding$sample.coords
  feature.coords <- swne.embedding$feature.coords

  feature.scores <- as.numeric(feature.scores[rownames(sample.coords)])
  feature.quantiles <- quantile(feature.scores, probs = quantiles)

  min.quantile <- feature.quantiles[[1]]; max.quantile <- feature.quantiles[[2]];
  feature.scores[feature.scores < min.quantile] <- min.quantile
  feature.scores[feature.scores > max.quantile] <- max.quantile
  feature.scores <- normalize_vector(feature.scores, method = "bounded")

  sample.coords$pt.size <- pt.size
  sample.coords$feature <- feature.scores

  if (!is.null(samples.plot)) {
    sample.coords <- sample.coords[samples.plot,]
    sample.coords <- droplevels(sample.coords)
  }

  ## Plot sample coordinates
  ggobj <- ggplot() +
    geom_point(data = sample.coords, aes(x, y, colour = feature),
               alpha = alpha.plot, size = pt.size) +
    theme_classic() + theme(axis.title = element_blank(), axis.ticks = element_blank(),
                            axis.line = element_blank(), axis.text = element_blank()) +
    scale_color_distiller(palette = color.palette, direction = 1, guide =
                            guide_colorbar(title = feature.name, ticks = F, label = F))

  if (!is.null(feature.coords)) {
    ggobj <- ggobj + geom_point(data = feature.coords, aes(x, y), size = 2.5, color = "darkred")
  }

  ## Plot factors
  if (nrow(H.coords.plot) > 0) {
    label.pts <- rbind(H.coords.plot, feature.coords)
    ggobj <- ggobj + geom_point(data = H.coords.plot, aes(x, y), size = 2.5, color = "darkblue")
    ggobj <- ggobj + ggrepel::geom_text_repel(data = label.pts, mapping = aes(x, y, label = name),
                                              size = label.size, box.padding = 0.15)
  } else if (!is.null(feature.coords)) {
    ggobj <- ggobj + ggrepel::geom_text_repel(data = feature.coords, mapping = aes(x, y, label = name),
                                              size = label.size, box.padding = 0.15)
  }

  return(ggobj)
}


#' Plots 2D embedding
#'
#' @param dim.scores 2D embedding coordinates. Must be N x 2 samples
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
#' @import ggrepel
#' @export
#'
PlotDims <- function(dim.scores, sample.groups = NULL, x.lab = "tsne1", y.lab = "tsne2",
                     main.title = NULL, pt.size = 1.0, font.size = 12, alpha.plot = 1.0, do.label = T,
                     label.size = 4, show.legend = T, show.axes = T, seed = NULL) {

  set.seed(seed)
  if (!is.null(sample.groups)) {
    sample.groups <- factor(sample.groups, levels = sample(levels(sample.groups)))
  }

  gg.df <- data.frame(x = dim.scores[,1], y = dim.scores[,2], sample.groups = sample.groups)
  ggobj <- ggplot(gg.df, aes(x, y)) + geom_point(size = pt.size, alpha = alpha.plot, aes(colour = sample.groups)) +
    theme(text = element_text(size = font.size)) +
    xlab(x.lab) + ylab(y.lab) + ggtitle(main.title) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = label.size)))

  if (!show.axes) {
    ggobj <- ggobj + theme_void()
  } else {
    ggobj <- ggobj + theme_classic()
  }

  if (do.label && is.factor(sample.groups)) {
    group.pts.x <- tapply(gg.df$x, gg.df$sample.groups, function(v) median(v))
    group.pts.y <- tapply(gg.df$y, gg.df$sample.groups, function(v) median(v))
    group.pts <- data.frame(x = group.pts.x, y = group.pts.y)
    group.pts$ident <- levels(gg.df$sample.groups)

    ggobj <- ggobj + ggrepel::geom_text_repel(data = group.pts, mapping = aes(label = ident), size = label.size,
                                              box.padding = 0.15)
  }

  if (!show.legend) { ggobj <- ggobj + theme(legend.position = "none") }

  ggobj
}


#' Plots 2d embedding with feature overlayed
#'
#' @param dim.scores 2D embedding coordinates. Must be N x 2 samples
#' @param feature.score Feature vector to overlay
#' @param feature.name Name of feature
#' @param x.lab X axis label
#' @param y.lab Y axis label
#' @param alpha.plot Data point transparency
#' @param quantiles Quantiles to trim outliers from
#' @param show.axes Plot x and y axes
#' @param pt.size Sample point size
#' @param font.size Font size for axis labels
#' @param color.palette RColorbrewer palette to use
#'
#' @return ggplot2 object with dim plot with feature overlayed
#'
#' @export
#'
FeaturePlotDims <- function(dim.scores, feature.scores, feature.name = NULL, x.lab = "tsne1", y.lab = "tsne2",
                            alpha.plot = 0.5, quantiles = c(0.01, 0.99), show.axes = T, pt.size = 1,
                            font.size = 12, color.palette = "YlOrRd") {

  sample.coords <- data.frame(x = dim.scores[,1], y = dim.scores[,2])

  feature.scores <- as.numeric(feature.scores[rownames(dim.scores)])
  feature.quantiles <- quantile(feature.scores, probs = quantiles)

  min.quantile <- feature.quantiles[[1]]; max.quantile <- feature.quantiles[[2]];
  feature.scores[feature.scores < min.quantile] <- min.quantile
  feature.scores[feature.scores > max.quantile] <- max.quantile
  feature.scores <- normalize_vector(feature.scores, method = "bounded")

  sample.coords$pt.size <- pt.size
  sample.coords$feature <- feature.scores

  ## Plot sample coordinates
  ggobj <- ggplot() +
    geom_point(data = sample.coords, aes(x, y, colour = feature),
               alpha = alpha.plot, size = pt.size) +
    xlab(x.lab) + ylab(y.lab) +
    scale_color_distiller(palette = color.palette, direction = 1, guide =
                            guide_colorbar(title = feature.name, ticks = F, label = F))

  if (!show.axes) {
    ggobj <- ggobj + theme_void()
  } else {
    ggobj <- ggobj + theme_classic()
  }

  ggobj <- ggobj + theme(text = element_text(size = font.size))

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
#' @importFrom reshape melt
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
  melt.m = cbind(rowInd = rep(1:rows, times = cols), colInd = rep(1:cols, each = rows), melt(m))
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


#' Extracts the exact colors used to plot each cluster (the hex codes) for a given color seed
#'
#' @param swne.embedding SWNE embedding (list of factor and sample coordinates) from EmbedSWNE
#' @param sample.groups Factor defining sample groups
#' @param seed Seed for sample groups color reproducibility
#'
#' @return vector of color hex codes with clusters as the name
#'
#' @export
ExtractSWNEColors <- function(swne.embedding, sample.groups, seed) {
  swne.ggobj <- PlotSWNE(swne.embedding, sample.groups = sample.groups, seed = seed)
  swne.ggobj <- ggplot2::ggplot_build(swne.ggobj)

  sample.colors <- swne.ggobj$data[[1]]$fill
  color.mapping <- unique(sample.colors)
  names(color.mapping) <- unique(as.character(sample.groups))

  return(color.mapping)
}
