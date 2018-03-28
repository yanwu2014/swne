#' Plots swne embedding with RNA velocity arrows overlaid
#'
#' @param swne.embedding SWNE embedding (list of factor and sample coordinates) from EmbedSWNE
#' @param swne.arrows DataFrame of Veloctyo arrow info from VeloctyoArrows
#' @param alpha.plot Data point transparency
#' @param sample.groups Factor defining sample groups
#' @param arrow.lwd Arrow width
#' @param arrow.alpha Arrow transparency
#' @param head.size Arrowhead size
#' @param do.label Label the sample groups
#' @param label.size Label font size
#' @param pt.size Sample point size
#' @param samples.plot Vector of samples to plot. Default is NULL, which plots all samples.
#' @param show.legend If sample groups defined, show legend
#' @param seed Seed for sample groups color reproducibility
#'
#' @return ggplot2 object with swne plot and RNA velocity arrows
#'
#' @export
#'
PlotSWNEVelocyto <- function(swne.embedding, swne.arrows, alpha.plot = 0.25, sample.groups = NULL,
                             arrow.lwd = 0.5, arrow.alpha = 0.5, head.size = 6e-3, do.label = F,
                             label.size = 4.5, pt.size = 1, samples.plot = NULL, show.legend = T,
                             seed = NULL) {
  ## Create SWNE plot
  swne.ggobj <- PlotSWNE(swne.embedding, alpha.plot, sample.groups, do.label = do.label,
                         label.size = label.size, pt.size = pt.size, samples.plot = samples.plot,
                         show.legend = show.legend, seed = seed)

  ## Add pre-computed arrows
  swne.ggobj <- swne.ggobj +
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, alpha = arrow.alpha), size = arrow.lwd,
                 arrow = grid::arrow(angle = 50, length = grid::unit(head.size, "npc")),
                 data = swne.arrows)
  return(swne.ggobj)
}


#' Function for computing arrows to add onto a ggplot2 object.
#' Adapted from velocyto.R
#'
#' @param emb 2D embedding (SWNE, t-SNE, etc.) coordinates. Must be N x 2 samples
#' @param vel RNA velocity object
#' @param pca.red PCA reduction. Must be N x nPCs samples.
#'
#' @export
#'
VelocytoArrows <- function(emb, vel, pca.red = NULL, n = 200, scale = "sqrt", corr.sigma = 0.05,
                           show.grid.flow = T, grid.n = 30, grid.sd = NULL, min.grid.cell.mass = 0.5,
                           min.arrow.size = NULL, arrow.scale = 1, max.grid.arrow.length = NULL,
                           n.cores = 16, diffusion.steps = 10, cc = NULL) {
  if (!requireNamespace("velocyto.R", quietly = TRUE)) {
    stop("Package \"velocyto.R\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  em <- as.matrix(vel$current);
  ccells <- intersect(rownames(emb),colnames(em));
  em <- em[,ccells]; emb <- emb[ccells,]
  nd <- as.matrix(vel$deltaE[,ccells])
  cgenes <- intersect(rownames(em),rownames(nd));
  nd <- nd[cgenes,]; em <- em[cgenes,]
  #vg <- rownames(em) %in% rownames(r)

  if(is.null(cc)) {
    # cosine projections
    cat("delta projections ... ")

    if(scale=='log') {
      cat("log ")
      cc <- velocyto.R:::colDeltaCorLog10(em,(log10(abs(nd)+1)*sign(nd)),nthreads=n.cores);
    } else if(scale=='sqrt') {
      cat("sqrt ")
      cc <- velocyto.R:::colDeltaCor(em,(sqrt(abs(nd))*sign(nd)),nthreads=n.cores);
    } else if(scale=='rank') {
      cat("rank ")
      cc <- velocyto.R:::colDeltaCor((apply(em,2,rank)),(apply(nd,2,rank)),nthreads=n.cores);
    } else { # linear
      cat("linear ")
      cc <- velocyto.R:::colDeltaCor(em,nd,nthreads=n.cores);
    }
    colnames(cc) <- rownames(cc) <- colnames(em)
    diag(cc) <- 0;
  }

  cat("knn ... ")
  if(n>nrow(cc)) { n <- nrow(cc) }
  # TODO: add kNN based on high-dimensional correlation or Euclidean distances
  # define kNNs based on the embedding (L2 distance)
  if (is.null(pca.red)) {
    emb.knn <- velocyto.R:::balancedKNN(t(emb),k=n,maxl=nrow(emb),dist='euclidean',n.threads=n.cores)
  } else {
    emb.knn <- velocyto.R:::balancedKNN(t(pca.red),k=n,maxl=nrow(emb),dist='euclidean',n.threads=n.cores)
  }

  diag(emb.knn) <- 1
  # caluclate transition probabilities (from col to row)
  cat("transition probs ... ")
  tp <- exp(cc/corr.sigma)*emb.knn
  #diag(tp) <- 0; #  should we allow the self-corelation for scaling?
  tp <- t(t(tp)/Matrix::colSums(tp))
  cat("done\n")

  if(show.grid.flow) {
    # show grid summary of the arrows
    # arrow estimates for each cell

    cat("calculating arrows ... ")
    arsd <- data.frame(t(velocyto.R:::embArrows(emb,as.matrix(tp),arrow.scale,n.cores)))
    ars <- data.frame(cbind(emb,emb+arsd));
    colnames(ars) <- c('x0','y0','x1','y1')
    colnames(arsd) <- c('xd','yd')

    # set up a grid
    cat("grid estimates ... ")
    rx <- range(c(range(ars$x0),range(ars$x1)))
    ry <- range(c(range(ars$y0),range(ars$y1)))
    gx <- seq(rx[1],rx[2],length.out=grid.n)
    gy <- seq(ry[1],ry[2],length.out=grid.n)

    # for each grid point calculate Gaussian-weighted delta average
    if(is.null(grid.sd)) {
      grid.sd <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)/2
      cat("grid.sd=",grid.sd," ")
    }
    if(is.null(min.arrow.size)) {
      min.arrow.size <- sqrt((gx[2]-gx[1])^2 + (gy[2]-gy[1])^2)*1e-2;
      cat("min.arrow.size=",min.arrow.size," ")
    }
    if(is.null(max.grid.arrow.length)) {
      max.grid.arrow.length <- sqrt(sum((par('pin')/c(length(gx),length(gy)))^2))*0.25
      cat("max.grid.arrow.length=",max.grid.arrow.length," ")
    }

    emb.arrows <- NULL
    for (x in gx) {
      # cell distances (rows-cells,columsn - grid points)
      cd <- sqrt(outer(emb[,2],-gy,'+')^2 + (x-emb[,1])^2)
      cw <- dnorm(cd, sd = grid.sd)
      # calculate x and y delta expectations
      gw <- Matrix::colSums(cw)
      cws <- pmax(1,Matrix::colSums(cw));
      gxd <- Matrix::colSums(cw*arsd$xd)/cws
      gyd <- Matrix::colSums(cw*arsd$yd)/cws

      al <- sqrt(gxd^2+gyd^2);
      vg <- gw >= min.grid.cell.mass & al >= min.arrow.size

      if(any(vg)) {
        col.arrows.df <- data.frame(x1 = rep(x, sum(vg)), y1 = gy[vg], x2 = x + gxd[vg], y2 = gy[vg] + gyd[vg])
        emb.arrows <- rbind(emb.arrows, col.arrows.df)
      }
    }
    cat("done\n")

  } else {
    ## calculate arrows, draw
    emb.arrows <- data.frame(x1 = numeric(nrow(emb)), y1 = numeric(nrow(emb)),
                             x2 = numeric(nrow(emb)), y2 = numeric(nrow(emb)))
    for(i in 1:nrow(emb)) {
      ## normalized directions to each point
      di <- t(t(emb)-emb[i,])
      di <- di/sqrt(Matrix::rowSums(di^2))*arrow.scale; di[i,] <- 0;
      di <- Matrix::colSums(di*tp[,i]) - Matrix::colSums(di*(tp[,i]>0)/sum(tp[,i]>0)); # relative to expected kNN center

      emb.arrows[i,"x1"] <- emb[colnames(em)[i],1]
      emb.arrows[i,"y1"] <- emb[colnames(em)[i],2]
      emb.arrows[i,"x2"] <- emb[colnames(em)[i],1] + di[1]
      emb.arrows[i,"y2"] <- emb[colnames(em)[i],2] + di[2]
    }
    cat("done\n")
  }

  return(emb.arrows)
}
