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
#' @importFrom proxy dist
#' @import umap
#'
get_factor_coords <- function(H, method, distance, n.neighbors, min.dist) {
  H <- t(H)
  distance <- tolower(distance)
  if(distance == "cor" || distance == "correlation") distance <- "pearson"
  if(distance == "mutual" || distance == "information" || distance == "mutual information" || distance == "ic") distance <- "IC"
  if (!distance %in% c("pearson", "IC", "cosine", "euclidean")) {
    stop(paste(c("Distance must be one of:", paste(c("pearson", "IC", "cosine", "euclidean"), collapse = ", ")), collapse = " "))
  }

  if (method == "sammon") {
    if (distance == "pearson") {
      H.dist <- sqrt(2*(1 - cor(H)))
    } else if (distance == "IC") {
      H.dist <- sqrt(2*(1 - as.matrix(usedist::dist_make(t(H), distance_fcn = MutualInf, method = "IC"))))
    } else if (distance == "cosine") {
      # H.dist <- sqrt(2*(1 - proxy::simil(H, method = "cosine", by_rows = F)))
      H.dist <- 1 - proxy::simil(H, method = "cosine", by_rows = F)
    } else if (distance == "euclidean") {
      H.dist <- proxy::dist(H, method = "Euclidean", by_rows = F)
    }
    H.coords <- MASS::sammon(H.dist, k = 2, niter = 250)$points
    rownames(H.coords) <- colnames(H)

  } else {
    stop("Invalid factor projection method")
  }

  H.coords <- apply(H.coords, 2, normalize_vector, method = "bounded")
  colnames(H.coords) <- c("x","y")

  return(H.coords)
}


#' Calculate sample coordinates using the NMF scores and NMF factor coordinates
#'
#' @import compiler
get_sample_coords <- function(H, H.coords, alpha, n_pull) {
  if (n_pull < 3) { n_pull <- 3 }
  if (is.null(n_pull) || n_pull > nrow(H.coords)) { n_pull <- nrow(H.coords) }

  sample.coords <- t(sapply(1:ncol(H), function(i) {
    pull.ix <- order(H[,i], decreasing = T)[1:n_pull]
    pull.sum <- sum(H[pull.ix,i]^alpha)
    x <- sum(H.coords[pull.ix,1] * ( H[pull.ix,i]^alpha/pull.sum ))
    y <- sum(H.coords[pull.ix,2] * ( H[pull.ix,i]^alpha/pull.sum ))
    return(c(x,y))
  }))
  colnames(sample.coords) <- c("x","y")
  rownames(sample.coords) <- colnames(H)
  return(sample.coords)
}; get_sample_coords <- cmpfun(get_sample_coords);



#' Embeds NMF factors and samples in 2 dimensions
#'
#' @param H NMF factors (factors x samples)
#' @param SNN Shared nearest neighbors matrix (or other similarity matrix)
#' @param alpha.exp Increasing alpha.exp increases how much the NMF factors "pull" the samples
#' @param snn.exp Decreasing snn.exp increases the effect of the similarity matrix on the embedding
#' @param n_pull Number of factors pulling on each sample. Must be >= 3
#' @param proj.method Method to use for projecting factors. Currently only supports "sammon"
#' @param dist.use Similarity function to use for calculating factor positions. Options include pearson (correlation), IC (mutual information), cosine, euclidean.
#'
#' @return A list of factor (H.coords) and sample coordinates (sample.coords) in 2D
#'
#' @export
#'
EmbedSWNE <- function(H, SNN = NULL, alpha.exp = 1, snn.exp = 1.0, n_pull = NULL,
                      proj.method = "sammon", dist.use = "cosine",
                      snn.factor.proj = T) {

  if (!is.null(SNN)) {
    if (!(nrow(SNN) == ncol(H) && ncol(SNN) == ncol(H))) {
      stop("Check the dimensions of your SNN. nrow(H) must equal nrow(SNN) and ncol(SNN)")
    }
  }

  if (nrow(H) < 4 && proj.method == "umap") {
    warning("Need at least 4 factors to use UMAP projection, defaulting to Sammon mapping")
    proj.method = "sammon"
  }

  ix <- colSums(H) > 0
  H <- H[,ix]

  if (!is.null(SNN)) {
    SNN <- SNN[ix, ix]
    SNN@x <- SNN@x^snn.exp
    SNN <- SNN/Matrix::rowSums(SNN)
  }

  if (!is.null(SNN) && snn.factor.proj) {
    H.smooth <- t(as.matrix(SNN %*% t(H)))
    H.coords <- get_factor_coords(H.smooth, method = proj.method,
                                  distance = dist.use)
  } else {
    H.coords <- get_factor_coords(H, method = proj.method,
                                  distance = dist.use)
  }

  H.coords <- data.frame(H.coords)
  H.coords$name <- rownames(H.coords); rownames(H.coords) <- NULL;

  sample.coords <- get_sample_coords(H, H.coords, alpha = alpha.exp, n_pull = n_pull)
  if (!is.null(SNN)) {
    sample.coords <- as.matrix(SNN %*% sample.coords)
    rownames(sample.coords) <- colnames(H)
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
#' @param overwrite Whether or not to overwrite any existing feature embedding
#'
#' @return swne.embedding with feature coordinates (feature.coords)
#'
#' @export
#'
EmbedFeatures <- function(swne.embedding, feature.assoc, features.embed, alpha.exp = 1, n_pull = NULL,
                          scale.cols = T, overwrite = T) {
  feature.assoc <- t(feature.assoc[features.embed,])
  stopifnot(nrow(swne.embedding$H.coords) == nrow(feature.assoc))

  if (scale.cols) {
    feature.assoc <- apply(feature.assoc, 2, normalize_vector, method = "bounded")
  }

  feature.coords <- get_sample_coords(feature.assoc, swne.embedding$H.coords, alpha = alpha.exp, n_pull = n_pull)
  feature.coords <- data.frame(feature.coords)
  feature.coords$name <- rownames(feature.coords); rownames(feature.coords) <- NULL;

  if (overwrite || is.null(swne.embedding$feature.coords)) {
    swne.embedding$feature.coords <- feature.coords
  } else {
    swne.embedding$feature.coords <- rbind(swne.embedding$feature.coords, feature.coords)
  }

  return(swne.embedding)
}


#' Projects new data onto an existing swne embedding
#'
#' @param swne.embedding Existing swne embedding from EmbedSWNE
#' @param H.test Test factor scores
#' @param SNN Test SNN matrix
#' @param alpha.exp Increasing alpha.exp increases how much the factors "pull" the samples
#' @param snn.exp Decreasing snn.exp increases the effect of the similarity matrix on the embedding
#' @param n_pull Number of factors pulling on each sample. Must be >= 3
#'
#' @return SWNE embedding for the test data
#'
#' @export
#'
ProjectSWNE <- function(swne.embedding, H.test, SNN = NULL, alpha.exp = 1, snn.exp = 1, n_pull = NULL) {
  sample.coords.test <- get_sample_coords(H.test, swne.embedding$H.coords, alpha = alpha.exp,
                                          n_pull = n_pull)

  if (!is.null(SNN)) {
    SNN <- SNN[rownames(sample.coords.test), rownames(swne.embedding$sample.coords)]

    snn.diag <- Matrix::Diagonal(ncol(H.test))
    rownames(snn.diag) <- colnames(snn.diag) <- colnames(H.test)
    SNN <- cbind(snn.diag, SNN); rm(snn.diag);

    SNN <- SNN^snn.exp
    SNN <- SNN/Matrix::rowSums(SNN)

    sample.coords.train <- as.matrix(swne.embedding$sample.coords)
    sample.coords.all <- rbind(sample.coords.test, sample.coords.train)

    sample.coords.test <- as.matrix(SNN %*% sample.coords.all)
    rownames(sample.coords.test) <- colnames(H.test)
    sample.coords.test <- data.frame(sample.coords.test)
  } else {
    sample.coords.test <- data.frame(sample.coords.test)
  }

  swne.embedding.test <- swne.embedding
  swne.embedding.test$sample.coords <- sample.coords.test
  return(swne.embedding.test)
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
  if (set.empty) {
    new.names[new.names == old.names] <- ""
  }
  swne.embedding$H.coords$name <- new.names
  return(swne.embedding)
}



#' Plots swne embedding
#'
#' @param swne.embedding SWNE embedding (list of factor and sample coordinates) from EmbedSWNE
#' @param sample.groups Factor defining sample groups
#' @param alpha.plot Data point transparency
#' @param do.label Label the sample groups
#' @param label.size Label font size
#' @param pt.size Sample point size
#' @param samples.plot Vector of samples to plot. Default is NULL, which plots all samples
#' @param show.legend If sample groups defined, show legend
#' @param seed Seed for sample groups color reproducibility
#' @param colors.use Vector of hex colors for each sample group. Vector names must align with sample.groups
#' @param use.brewer.pal Use RColorBrewer 'Paired' palette instead of default ggplot2 color wheel
#' @param contour.geom Plot contour as either a line or a filled in region
#' @param contour.alpha Transparency of the contour line/region
#'
#' @return ggplot2 object with swne plot
#'
#' @import ggrepel
#' @import RColorBrewer
#' @export
#'
PlotSWNE <- function(swne.embedding, sample.groups, alpha.plot = 0.25, do.label = F,
                     label.size = 4.5, pt.size = 1, samples.plot = NULL, show.legend = T, seed = NULL,
                     colors.use = NULL, use.brewer.pal = F, contour.geom = "path",
                     contour.alpha = 0.25) {
  H.coords <- swne.embedding$H.coords
  H.coords.plot <- subset(H.coords, name != "")

  sample.coords <- swne.embedding$sample.coords
  feature.coords <- swne.embedding$feature.coords
  sample.coords$pt.size <- pt.size

  sample.groups <- factor(sample.groups)
  if (!is.null(names(sample.groups))) {
    samples.plot <- intersect(names(sample.groups), rownames(sample.coords))
    sample.groups <- sample.groups[samples.plot]
    sample.coords <- sample.coords[samples.plot,]
  } else {
    if (length(sample.groups) != nrow(sample.coords)) {
      stop("length(sample.groups) must be equal to nrow(swne.embedding$sample.coords)")
    }
  }

  n.missing <- sum(is.na(sample.groups))
  if (n.missing > 0) {
    na.ix <- is.na(sample.groups)
    label.ix <- !is.na(sample.groups)
    na.sample.coords <- sample.coords[na.ix,]
    sample.coords <- sample.coords[label.ix,]
    sample.groups <- sample.groups[label.ix]
  }
  set.seed(seed)
  sample.groups <- factor(sample.groups, levels = sample(levels(sample.groups)))
  sample.coords$sample.groups <- sample.groups

  if (!is.null(samples.plot)) {
    sample.coords <- sample.coords[samples.plot,]
    sample.coords <- droplevels(sample.coords)
  }

  ## Create ggplot2 object
  ggobj <- ggplot()

  ## Plot confidence ellipse
  if (!is.null(swne.embedding$contour.data)) {
    if(!contour.geom %in% c("path", "polygon")) stop("contour.geom must be one of: 'path', 'polygon'")
    ggobj <- ggobj +
      stat_ellipse(mapping = aes(x, y), data = swne.embedding$contour.data,
                   geom = contour.geom, position = "identity", type = "t", level = 0.95,
                   na.rm = T, inherit.aes = F, linetype = 2, alpha = contour.alpha,
                   colour = "darkred")
  }

  if (n.missing > 0) {
    ## Plot samples with missing labels
    ggobj <- ggobj +
      geom_point(data = na.sample.coords, aes(x, y), alpha = alpha.plot, size = pt.size, color = "lightgrey")
  }

  ## Plot sample coordinates
  ggobj <- ggobj +
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
    label.pts <- subset(label.pts, name != "NA")

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
  } else {
    ggobj <- ggobj + theme(legend.title = element_blank())
  }

  if (use.brewer.pal) {
    brewer.colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(nlevels(sample.groups))
    names(brewer.colors) <- levels(sample.groups)
    ggobj <- ggobj + scale_color_manual(values = brewer.colors)
  }

  if (!is.null(colors.use)) {
    if(!all(levels(sample.groups) %in% names(colors.use))) {
      stop("Must specify colors for each group")
    }
    ggobj <- ggobj + scale_color_manual(values = colors.use)
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
#' @param contour.geom Plot contour as either a line or a filled in region
#' @param contour.alpha Transparency of the contour line/region
#'
#' @return ggplot2 object with swne plot with feature overlayed
#'
#' @import ggrepel
#' @export
#'
FeaturePlotSWNE <- function(swne.embedding, feature.scores, feature.name = NULL, alpha.plot = 0.5,
                            quantiles = c(0.01, 0.99), samples.plot = NULL, label.size = 4.5,
                            pt.size = 1, color.palette = "YlOrRd", contour.geom = "path",
                            contour.alpha = 0.25) {
  H.coords <- swne.embedding$H.coords
  H.coords.plot <- subset(H.coords, name != "")
  sample.coords <- swne.embedding$sample.coords
  feature.coords <- swne.embedding$feature.coords

  feature.scores <- feature.scores[!is.na(feature.scores)]
  sample.coords <- sample.coords[names(feature.scores),]
  feature.scores <- as.numeric(feature.scores[rownames(sample.coords)])
  feature.quantiles <- quantile(feature.scores, probs = quantiles)

  min.quantile <- feature.quantiles[[1]]; max.quantile <- feature.quantiles[[2]];
  feature.scores[feature.scores < min.quantile] <- min.quantile
  feature.scores[feature.scores > max.quantile] <- max.quantile
  feature.scores <- normalize_vector(feature.scores, method = "bounded")

  sample.coords$pt.size <- pt.size
  sample.coords$feature <- feature.scores
  sample.coords <- sample.coords[order(sample.coords$feature),]

  if (!is.null(samples.plot)) {
    sample.coords <- sample.coords[samples.plot,]
    sample.coords <- droplevels(sample.coords)
  }

  ## Plot sample coordinates
  ggobj <- ggplot() +
    geom_point(data = sample.coords, aes(x, y, colour = feature),
               alpha = alpha.plot, size = pt.size) +
    theme_classic() + theme(axis.title = element_blank(), axis.ticks = element_blank(),
                            axis.line = element_blank(), axis.text = element_blank(),
                            legend.text = element_text(size = 12)) +
    scale_color_distiller(palette = color.palette, direction = 1, guide =
                            guide_colorbar(title = feature.name, ticks = T, label = T))

  if (!is.null(feature.coords)) {
    ggobj <- ggobj + geom_point(data = feature.coords, aes(x, y), size = 2.5, color = "darkred")
  }

  ## Plot confidence ellipse
  if (!is.null(swne.embedding$contour.data)) {
    if(!contour.geom %in% c("path", "polygon")) stop("contour.geom must be one of: 'path', 'polygon'")
    ggobj <- ggobj +
      stat_ellipse(mapping = aes(x, y), data = swne.embedding$contour.data,
                   geom = contour.geom, position = "identity", type = "t", level = 0.95,
                   na.rm = T, inherit.aes = F, linetype = 2, alpha = contour.alpha,
                   colour = "darkred")
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
#' @param colors.use Vector of hex colors for each sample group. Vector names must align with sample.groups
#' @param use.brewer.pal Use RColorBrewer 'Paired' palette instead of default ggplot2 color wheel
#'
#' @return ggplot2 object with 2d plot
#'
#' @import ggrepel
#' @import RColorBrewer
#' @export
#'
PlotDims <- function(dim.scores, sample.groups, x.lab = "tsne1", y.lab = "tsne2",
                     main.title = NULL, pt.size = 1.0, font.size = 12, alpha.plot = 1.0, do.label = T,
                     label.size = 4, show.legend = T, show.axes = T, seed = NULL, colors.use = NULL,
                     use.brewer.pal = F) {

  sample.groups <- factor(sample.groups)
  if (!is.null(names(sample.groups))) {
    samples.plot <- intersect(names(sample.groups), rownames(dim.scores))
    sample.groups <- sample.groups[samples.plot]
    dim.scores <- dim.scores[samples.plot,]
  } else {
    if (length(sample.groups) != nrow(dim.scores)) {
      stop("length(sample.groups) must be equal to nrow(dim.scores)")
    }
  }

  n.missing <- sum(is.na(sample.groups))
  if (n.missing > 0) {
    na.ix <- is.na(sample.groups)
    label.ix <- !is.na(sample.groups)
    na.dim.scores <- dim.scores[na.ix,]
    dim.scores <- dim.scores[label.ix,]
    sample.groups <- sample.groups[label.ix]
  }
  set.seed(seed)
  sample.groups <- factor(sample.groups, levels = sample(levels(sample.groups)))

  gg.df <- data.frame(x = dim.scores[,1], y = dim.scores[,2], sample.groups = sample.groups)

  ggobj <- ggplot()
  if (n.missing > 0) {
    ## Plot samples with missing labels
    na.gg.df <- data.frame(x = na.dim.scores[,1], y = na.dim.scores[,2])
    ggobj <- ggobj + geom_point(data = na.gg.df, aes(x, y), alpha = alpha.plot,
                                size = pt.size, color = "lightgrey")
  }

  ggobj <- ggobj + geom_point(data = gg.df, size = pt.size, alpha = alpha.plot,
                              aes(x, y, colour = sample.groups)) +
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
    group.pts <- subset(group.pts, ident != "NA")

    ggobj <- ggobj + ggrepel::geom_text_repel(data = group.pts, mapping = aes(x, y, label = ident),
                                              size = label.size, box.padding = 0.15)
  }

  if (!show.legend) {
    ggobj <- ggobj + theme(legend.position = "none")
  } else {
    ggobj <- ggobj + theme(legend.title = element_blank())
  }

  if (use.brewer.pal) {
    brewer.colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(nlevels(sample.groups))
    names(brewer.colors) <- levels(sample.groups)
    ggobj <- ggobj + scale_color_manual(values = brewer.colors)
  }

  if (!is.null(colors.use)) {
    if(!all(levels(sample.groups) %in% names(colors.use))) {
      stop("Must specify colors for each group")
    }
    ggobj <- ggobj + scale_color_manual(values = colors.use)
  }

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
#' @import ggplot2
#' @import ggrepel
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
  sample.coords <- sample.coords[order(sample.coords$feature),]


  ## Plot sample coordinates
  ggobj <- ggplot() +
    geom_point(data = sample.coords, aes(x, y, colour = feature),
               alpha = alpha.plot, size = pt.size) +
    xlab(x.lab) + ylab(y.lab) +
    scale_color_distiller(palette = color.palette, direction = 1, guide =
                            guide_colorbar(title = feature.name, ticks = F, label = T))

  if (!show.axes) {
    ggobj <- ggobj + theme_void()
  } else {
    ggobj <- ggobj + theme_classic()
  }

  ggobj <- ggobj + theme(text = element_text(size = font.size), legend.text = element_text(size = 12))

  ggobj
}


#' Plots a heatmap using ggplot2
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
#' @param dot.highlight.cutoff Heatmap cells above this value  will be highlighted with a dot
#'
#' @return ggplot2 heatmap
#'
#' @importFrom reshape melt
#'
#' @export
#'
ggHeat <- function(m, rescaling = 'none', clustering = 'none',
                   labCol = T, labRow = T, border = F,
                   heatscale = c(low = 'skyblue', mid = 'white', high = 'tomato'),
                   legend.title = NULL, x.lab.size = 10,
                   y.lab.size = 10,
                   dot.highlight.cutoff = Inf) {
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

  melt.m = reshape2::melt(m)
  g2 <- ggplot(data = melt.m) + geom_tile(aes(x = Var2, y = Var1, fill = value))


  ## get rid of grey panel background and gridlines
  g2 = g2 + theme(panel.grid.minor = element_line(colour = NA), panel.grid.major = element_line(colour=NA),
                  panel.background = element_rect(fill = NA, colour = NA), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = x.lab.size),
                  axis.ticks = element_blank(), axis.text.y = element_text(size = y.lab.size, hjust = 1, vjust = 0.5),
                  axis.title = element_blank())

  ## finally add the fill colour ramp of your choice (default is blue to red)-- and return
  g2 <- g2 + scale_fill_gradient2(low = heatscale[1], mid = heatscale[2], high = heatscale[3], guide = guide_colorbar(title = legend.title))

  ## Add dots to highlight certain cells
  dot.melt.m <- subset(melt.m, abs(value) > dot.highlight.cutoff)
  if (nrow(dot.melt.m) > 0) {
    g2 <- g2 + geom_point(aes(x = colInd, y = rowInd), data = dot.melt.m)
  }

  ## Remove axis labels
  if(labCol == F) g2 = g2 + theme(axis.text.y = element_blank())
  if(labRow == F) g2 = g2 + theme(axis.text.x = element_blank())

  return(g2)
}


#' Make simple barplot from vector
#'
#' @param x Named numeric vector
#' @param y.lim Max y-axis
#' @param fill.color Bar color
#'
#' @return barplot
#' @import ggplot2
#' @export
#'
ggBarplot <- function(x, y.lim = NULL, fill.color = "lightgrey") {
  unique.names <- paste0(names(x), 1:length(x))
  barplot.df <- data.frame(Y = x, X = factor(unique.names, levels = unique.names),
                           color = fill.color)

  if(length(fill.color) == 1) {
    ggobj <- ggplot(data = barplot.df, aes(x = X, y = Y)) +
      geom_bar(position = "dodge", stat = "identity", width = 0.8, color = "black", fill = fill.color) +
      scale_x_discrete(labels = names(x))
  } else {
    ggobj <- ggplot(data = barplot.df, aes(x = X, y = Y, fill = color)) +
      geom_bar(position = "dodge", stat = "identity", width = 0.8, color = "black") +
      scale_x_discrete(labels = names(x))
  }


  ggobj <- ggobj + theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 14, angle = 90, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"))

  if(!is.null(y.lim)) {
    ggobj <- ggobj + coord_cartesian(ylim = y.lim)
  }
  ggobj
}



#' Validates gene embeddings by plotting cluster logFC vs gene factor loading logFC
#' Warns users if embedded genes fall below minimum logFC threshold
#'
#' @param W Gene loadings matrix
#' @param norm.counts Normalized gene expression matrix
#' @param genes.embed Genes to embed onto SWNE plot
#' @param sample.groups Sample groupings
#' @param n.bins Number of bins for the hexbin plot
#' @param min.cluster.logfc Minimum cluster logFC for an embedded gene
#' @param min.factor.logfc Minimum factor logFC for an embedded gene
#' @param font.size Axis label font size
#' @param label.size Gene label size
#' @param eps Pseudocount to add to fold-change calculations
#'
#' @return Dataframe with cluster logFC and factor logFC. Generates scatterplot as a side effect.
#'
#' @import ggrepel
#' @export
#'
CheckGeneEmbedding <- function(W, norm.counts, genes.embed, sample.groups, n.bins = 50,
                               min.cluster.logfc = 1.5, min.factor.logfc = 1.5,
                               font.size = 12, label.size = 5, eps = 1e-4) {
  if(!all(rownames(W) == rownames(norm.counts))) {
    stop("Rownames of W must match rownames of norm.counts")
  }
  W <- W/colSums(W)

  gene.factor.logfc <- log2(apply(W, 1, function(x) {
    max.i <- which.max(x)
    (x[[max.i]] + eps)/(mean(x[-1*max.i]) + eps)
  }))

  gene.cell.logfc <- log2(apply(norm.counts, 1, function(x) {
    cl.avg <- tapply(x, sample.groups, mean)
    cl.max <- names(cl.avg[which.max(cl.avg)])
    (mean(x[sample.groups == cl.max]) + eps)/(mean(x[sample.groups != cl.max]) + eps)
  }))

  gene.logfc.df <- data.frame(cluster = gene.cell.logfc, factor = gene.factor.logfc)
  gene.logfc.df$embedded <- factor(rownames(gene.logfc.df) %in% genes.embed)
  gene.logfc.df <- gene.logfc.df[order(gene.logfc.df$embedded),]

  gg.obj <- ggplot() +
    geom_hex(data = gene.logfc.df, mapping = aes(factor, cluster), bins = n.bins, alpha = 1) +
    theme_classic() + theme(text = element_text(size = font.size)) +
    xlab("Max factor logFC") + ylab("Max cluster logFC") +
    scale_fill_gradient(low = "skyblue", high = "tomato3")

  gene.logfc.label <- subset(gene.logfc.df, embedded == "TRUE")
  gene.logfc.label$name <- rownames(gene.logfc.label)

  warning.genes <- rownames(subset(gene.logfc.label, cluster < min.cluster.logfc | factor < min.factor.logfc))
  if(length(warning.genes) > 0) {
    print("Warning, the following genes may not be good candidates for embedding:")
    print(warning.genes)
  }

  gg.obj <- gg.obj +
    geom_point(data = gene.logfc.label, mapping = aes(factor, cluster), size = 2, alpha = 1, color = "darkred") +
    ggrepel::geom_text_repel(data = gene.logfc.label, mapping = aes(factor, cluster, label = name),
                             size = label.size, box.padding = 0.15) +
    geom_vline(xintercept = min.factor.logfc, linetype = "dotted", size = 1.25, color = "black") +
    geom_hline(yintercept = min.cluster.logfc, linetype = "dotted", size = 1.25, color = "black")
  print(gg.obj)

  return(gene.logfc.df)
}



#' Plot decrease in reconstruction error versus random noise.
#' The point where the change in reconstruction error matches the change in reconstruction
#' error for the randomized matrix is the optimal number of factors to use.
#'
#' @param k.err.diff Reduction in reconstruction error above noise output by FindNumFactors
#' @param font.size Font size to use for plotting
#'
#' @import ggplot2
#' @export
#'
PlotFactorSelection <- function(k.err.diff, font.size = 12) {
  if(is.list(k.err.diff)) {
    k.err.diff <- k.err.diff[["err"]]
  }
  err.del.df <- data.frame(y = k.err.diff, x = factor(names(k.err.diff), levels = names(k.err.diff)))

  ggplot(data = err.del.df, aes(x, y)) +
    geom_line(aes(group = 1)) + geom_point(size = 2.0) +
    theme_classic() + xlab("Number of factors") + ylab("Error reduction above noise") +
    theme(axis.text.x = element_text(hjust = 1, size = font.size, angle = 90, color = "black"),
          axis.text = element_text(size = font.size, color = "black"),
          axis.title = element_text(size = font.size, color = "black"),
          legend.position = c(0.8,0.8),
          legend.text = element_text(size = font.size),
          legend.title = element_text(size = font.size)) +
    geom_hline(yintercept = 0.0, linetype = "dashed", color = "darkred")
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
