#### NMF initialization functions

## NNSVD helper functions
.pos <- function(x) { as.numeric(x >= 0) * x }
.neg <- function(x) { - as.numeric(x < 0) * x }
.norm <- function(x) { sqrt(drop(crossprod(x))) }


#' Nonnegative SVD initialization
#'
#' @import compiler
nnsvd_init <- function(A, k, LINPACK) {
  size <- dim(A);
  m <- size[1]; n <- size[2];

  W <- matrix(0, m, k);
  H <- matrix(0, k, n);

  #1st SVD --> partial SVD rank-k to the input matrix A.
  s = svd(A, k, k, LINPACK = LINPACK);
  U <- s$u; S <- s$d; V <- s$v

  #choose the first singular triplet to be nonnegative
  W[,1] = sqrt(S[1]) * abs(U[,1]);
  H[1,] = sqrt(S[1]) * abs(t(V[,1]));

  # second SVD for the other factors (see table 1 in Boutsidis' paper)
  for( i in seq(2,k) ) {
    uu = U[,i]; vv = V[,i];
    uup = .pos(uu); uun = .neg(uu) ;
    vvp = .pos(vv); vvn = .neg(vv);
    n_uup = .norm(uup);
    n_vvp = .norm(vvp) ;
    n_uun = .norm(uun) ;
    n_vvn = .norm(vvn) ;
    termp = n_uup %*% n_vvp; termn = n_uun %*% n_vvn;
    if (termp >= termn) {
      W[,i] = sqrt(S[i] * termp) * uup / n_uup;
      H[i,] = sqrt(S[i] * termp) * vvp / n_vvp;
    } else {
      W[,i] = sqrt(S[i] * termn) * uun / n_uun;
      H[i,] = sqrt(S[i] * termn) * vvn / n_vvn;
    }
  }

  return(list(W = W, H = H))
}; nnsvd_init <- compiler::cmpfun(nnsvd_init);


#' Independent component analysis initialization.
#' Negative values are set to a small, random number.
#'
#' @importFrom ica icafast
#'
ica_init <- function(A, k, ica.fast = F) {
  if (ica.fast) {
    pc.res.h <- irlba::irlba(t(A), nv = 100, maxit = 250, center = T)
    ica.res.h <- ica::icafast(pc.res.h$u, nc = k, maxit = 25, tol = 1e-4)
    return(list(W = (A - Matrix::rowMeans(A)) %*% ica.res.h$S,
                H = t(ica.res.h$S)))
  } else {
    ica.res <- ica::icafast(t(A), nc = k, maxit = 25, tol = 1e-4)
    return(list(W = ica.res$M, H = t(ica.res$S)))
  }
}


#' KL divergence. Pseudocounts added to avoid NAs
kl_div <- function(x, y, pseudocount = 1e-12) {
  x <- x + pseudocount; y <- y + pseudocount;
  stopifnot(length(x) == length(y))
  return(x*log(x/y) - x + y)
}



#' Runs NMF decomposition on a data matrix: A = WH.
#'
#' @param A Input data matrix
#' @param k Number of NMF factors
#' @param alpha Regularization parameter
#' @param init Initialization method: ica, nnsvd, or random.
#' @param n.cores Number of cores
#' @param loss Type of loss function to use
#' @param n.rand.init If random initialization is used, number of random restarts
#' @param init.zeros What to do with zeros in the initialization
#' @param max.iter Maximum number of iterations
#' @param ica.fast If using default ICA initialization, run PCA first to speed up ICA
#'
#' @return List of W (features x factors) and H (factors x samples)
#'
#' @import NNLM
#' @export
#'
RunNMF <- function(A, k, alpha = 0, init = "ica", n.cores = 1, loss = "mse",
                   init.zeros = "random", max.iter = 500, ica.fast = F) {
  if (any(A < 0)) stop('The input matrix contains negative elements !')
  if (k < 3) stop("k must be greater than or equal to 3 to create a viable SWNE plot")

  if (!init %in% c("ica", "nnsvd", "random")) {
    stop("Invalid initialization method")
  }

  if (!init.zeros %in% c("random", "uniform")) {
    stop("init.zeros should be either 'random', or 'uniform'")
  }

  A <- as.matrix(A)
  if (any(A < 0)) { stop("Input matrix has negative values") }

  if (init == "ica") {
    nmf.init <- ica_init(A, k, ica.fast = ica.fast)
  } else if (init == "nnsvd") {
    nmf.init <- nnsvd_init(A, k, LINPACK = T)
  } else {
    nmf.init <- NULL
  }

  if(is.null(nmf.init)) {
    # nmf.res.list <- lapply(1:n.rand.init, function(i) nnmf(A, k = k, alpha = alpha, init = nmf.init,
    #                                                        n.threads = n.cores, loss = loss,
    #                                                        max.iter = max.iter, verbose = F))
    # err <- sapply(nmf.res.list, function(x) tail(x[[loss]], n = 1))
    # nmf.res <- nmf.res.list[[which.min(err)]]
    nmf.res <- NNLM::nnmf(A, k = k, alpha = alpha, n.threads = n.cores, loss = loss, max.iter = max.iter)

  } else {
    ## Deal with zeros in the nmf initialization
    A.mean <- mean(A)
    zero.eps <- 1e-6

    nmf.init$W[nmf.init$W < zero.eps] <- 0; nmf.init$H[nmf.init$H < zero.eps] <- 0;
    zero.idx.w <- which(nmf.init$W == 0); zero.idx.h <- which(nmf.init$H == 0);

    if (init.zeros == "random") {
      nmf.init$W[zero.idx.w] <- runif(length(zero.idx.w), 0, A.mean/100)
      nmf.init$H[zero.idx.h] <- runif(length(zero.idx.h), 0, A.mean/100)
    } else if (init.zeros == "uniform") {
      nmf.init$W[zero.idx.w] <- A.mean/100
      nmf.init$H[zero.idx.h] <- A.mean/100
    }
    nmf.res <- NNLM::nnmf(A, k = k, alpha = alpha, init = nmf.init, n.threads = n.cores, loss = loss, max.iter = max.iter)
  }

  colnames(nmf.res$W) <- rownames(nmf.res$H) <- sapply(1:ncol(nmf.res$W), function(i) paste("factor", i, sep = "_"))
  return(nmf.res)
}


#' Determines the optimal number of NMF factors to use by comparing the reduction in
#' reconstruction error vs the reduction in reconstruction error for a randomized
#' matrix. The optimal number of factors is where the decrease in reconstruction
#' errors intersect
#'
#' Adapted from Frigyesi et al, 2008. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2623306/
#'
#' @param A Input data matrix
#' @param k.range Range of NMF factors to fit over
#' @param n.cores Number of threads
#' @param do.plot Whether to plot the reconstruction error
#' @param seed Random seed for randomizing matrix
#' @param loss Loss function to use for NMF
#' @param max.iter Maximum iterations for NMF run
#'
#' @return Reconstruction error at each number of NMF factors specified in k.range
#'
#' @import NNLM
#' @export
#'
FindNumFactors <- function(A, k.range = seq(2,12,2), n.cores = 1, do.plot = T,
                           seed = NULL, loss = "mse", max.iter = 250) {
  if (!loss %in% c("mse", "mkl")) { stop("Invalid loss function") }
  if (ncol(A) > 15000) print("Warning: This function can be slow for very large datasets")
  if (!is.null(seed)) { set.seed(seed) }
  if (length(k.range) < 3) { stop("k.range sequence must have at least 3 values") }

  A <- as.matrix(A)
  A.rand <- matrix(sample(A), nrow(A), ncol(A))
  k.err <- sapply(k.range, function(k) {
    z <- NNLM::nnmf(A, k, n.threads = n.cores, verbose = 0, max.iter = max.iter)
    z.rand <- NNLM::nnmf(A.rand, k, n.threads = n.cores, verbose = 0, max.iter = max.iter)

    A.hat <- with(z, W %*% H)
    A.hat.rand <- with(z.rand, W %*% H)

    if (loss == "mse") {
      err  <- mean((A.hat - A)^2)
      err.rand <- mean((A.hat.rand - A.rand)^2)
    } else if (loss == "mkl") {
      err <- mean(kl_div(A.hat, A))
      err.rand <- mean(kl_div(A.hat.rand, A.rand))
    }

    return(rbind(err, err.rand))
  })
  rownames(k.err) <- c("err", "err.rand")
  colnames(k.err) <- k.range


  if (do.plot) {
    print(PlotFactorSelection(k.err, font.size = 14))
  }

  min.idx <- which.max(k.err[match("err", rownames(k.err)), ])
  res <- list()
  res$err <- k.err
  res$k <- k.range[[min.idx]]

  return(res)
}



#' Projects new features onto existing factor decomposition
#'
#' @param newdata New data matrix
#' @param H Existing factor scores
#' @param alpha Regularization parameter
#' @param loss Loss function to use
#' @param n.cores Number of cores to use
#'
#' @return W matrix (feature loadings) for newdata
#'
#' @import NNLM
#' @export
#'
ProjectFeatures <- function(newdata, H, alpha = rep(0,3), loss = "mse", n.cores = 1) {
  lm.out <- NNLM::nnlm(t(H), t(newdata), alpha = alpha, loss = loss, n.threads = n.cores)
  return(t(lm.out$coefficients))
}


#' Projects new samples onto existing factor decomposition
#'
#' @param newdata New data matrix
#' @param W Existing feature loadings
#' @param features.use Subset of features to use (default is all features)
#' @param alpha Regularization parameter
#' @param loss Loss function to use
#' @param n.cores Number of cores to use
#'
#' @return H matrix (scores) for newdata
#'
#' @import NNLM
#' @export
#'
ProjectSamples <- function(newdata, W, features.use = NULL, alpha = rep(0,3), loss = "mse", n.cores = 1) {
  if (is.null(features.use)) {
    features.use <- intersect(rownames(newdata), rownames(W))
  }

  lm.out <- NNLM::nnlm(W[features.use,], newdata[features.use,], alpha = alpha, loss = loss, n.threads = n.cores)
  return(lm.out$coefficients)
}

