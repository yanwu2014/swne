#### NMF initialization functions

## NNSVD helper functions
.pos <- function(x) { as.numeric(x >= 0) * x }

.neg <- function(x) { - as.numeric(x < 0) * x }

.norm <- function(x) { sqrt(drop(crossprod(x))) }

#' Nonnegative SVD initialization
.nnsvd_init <- function(A, k, LINPACK, eps, init.zeros) {
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

  #actually these numbers are zeros
  W[W < eps] <- 0;
  H[H < eps] <- 0;

  ind1 <- W == 0; ind2 <- H == 0;
  A.mean <- mean(A);

  if (init.zeros == "random") {
    n1 <- sum(ind1); n2 <- sum(ind2);
    A.mean <- mean(A);
    W[ind1] <-  runif(n1, min = 0, max = A.mean) / 100
    H[ind2] <-  runif(n2, min = 0, max = A.mean) / 100
  } else if (init.zeros == "uniform") {
    W[ind1] <-  A.mean / 100
    H[ind2] <-  A.mean / 100
  }

  return(list(W = W, H = H))
}


#' Independent component analysis initialization.
#' Negative values are set to a small, random number.
#'
#' @import ica
#'
.ica_init <- function(A, k, eps, init.zeros) {
  ica.res <- ica::icafast(t(A), nc = k)
  nmf.init <- list(W = ica.res$M, H = t(ica.res$S))

  A.mean <- mean(A)
  nmf.init$W[nmf.init$W < eps] <- 0; nmf.init$H[nmf.init$H < eps] <- 0;
  zero.idx.w <- which(nmf.init$W == 0); zero.idx.h <- which(nmf.init$H == 0);

  if (init.zeros == "random") {
    nmf.init$W[zero.idx.w] <- runif(length(zero.idx.w), 0, A.mean/100)
    nmf.init$H[zero.idx.h] <- runif(length(zero.idx.h), 0, A.mean/100)
  } else if (init.zeros == "uniform") {
    nmf.init$W[zero.idx.w] <- A.mean / 100
    nmf.init$H[zero.idx.h] <- A.mean / 100
  }

  return(nmf.init)
}


#' KL divergence. Pseudocounts added to avoid NAs
.kl_div <- function(x, y, pseudocount = 1e-12) {
  x <- x + pseudocount; y <- y + pseudocount;
  stopifnot(length(x) == length(y))
  return(x*log(x/y) - x + y)
}


#' Determines the optimal number of NMF factors to use via reconstruction error
#'
#' @param A Input data matrix
#' @param k.range Range of NMF factors to fit over
#' @param alpha Regularization parameter
#' @param n.cores Number of threads
#' @param do.plot Whether to plot the reconstruction error
#' @param seed Random seed for selecting missing data
#' @param na.frac Fraction of data to set as missing
#' @param loss Loss function to use for NMF
#' @param recon.err Error function to minimize
#' @param max.iter Maximum iterations for NMF run
#'
#' @return Reconstruction error at each number of NMF factors specified in k.range
#'
#' @import NNLM
#' @export
#'
FindNumFactors <- function(A, k.range = seq(1,10,1), alpha = 0, n.cores = 1, do.plot = T,
                           seed = NULL, na.frac = 0.3, loss = "mse", recon.err = "mse", max.iter = 1000) {
  if (!is.null(seed)) { set.seed(seed) }
  if (!loss %in% c("mse", "mkl")) { stop("Invalid loss function") }
  if (!recon.err %in% c("mse", "mkl", "pearson", "spearman")) { stop("Invalid error function") }
  if (ncol(A) > 15000) print("Warning: This function can be slow for very large datasets")

  A <- as.matrix(A)
  nzero <- which(A > 0)
  # ind <- sample(nzero, na.frac*length(nzero));
  ind <- sample(length(A), na.frac*length(A))
  A2 <- A;
  A2[ind] <- NA;

  A.ind <- as.numeric(A[ind])
  err <- sapply(k.range, function(k) {
    z <- NNLM::nnmf(A2, k, alpha = c(alpha, alpha, 0), n.threads = n.cores, verbose = 0, loss = loss)
    A.hat <- with(z, W %*% H)
    A.hat.ind <- as.numeric(A.hat[ind])

    mse  <- mean((A.ind - A.hat.ind)^2)
    mkl <- mean(.kl_div(A.ind, A.hat.ind))
    pearson.r <- cor(A.ind, A.hat.ind)
    spearman.r <- cor(A.ind, A.hat.ind, method = "spearman")

    return(c(mse, mkl, pearson.r, spearman.r))
  })
  rownames(err) <- c("mse", "mkl", "pearson", "spearman")

  if (recon.err %in% c("pearson", "spearman")) {
    min.idx <- which.max(err[recon.err,])
  } else {
    min.idx <- which.min(err[recon.err,])
  }

  if (do.plot) {
    plot(k.range, err[recon.err,], col = "blue", type = 'b', main = "Model selection",
         xlab = "Number of factors", ylab = recon.err)
  }

  res <- list()
  res$err <- err
  res$k <- k.range[[min.idx]]

  return(res)
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
#'
#' @return List of W (features x factors) and H (factors x samples)
#'
#' @import NNLM
#' @export
#'
RunNMF <- function(A, k, alpha = 0, init = "random", n.cores = 1, loss = "mse", n.rand.init = 5,
                   init.zeros = "uniform") {
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
  A.mean <- mean(A)

  if (init == "ica") {
    nmf.init <- .ica_init(A, k, eps = 1e-8, init.zeros = init.zeros)
  } else if (init == "nnsvd") {
    nmf.init <- .nnsvd_init(A, k, LINPACK = T, eps = 1e-8, init.zeros = init.zeros)
  } else {
    nmf.init <- NULL
  }


  if(is.null(nmf.init)) {
    nmf.res.list <- lapply(1:n.rand.init, function(i) nnmf(A, k = k, alpha = alpha, init = nmf.init,
                                                           n.threads = n.cores, loss = loss, verbose = F))
    err <- sapply(nmf.res.list, function(x) tail(x[[loss]], n = 1))
    nmf.res <- nmf.res.list[[which.min(err)]]
  } else {
    nmf.res <- NNLM::nnmf(A, k = k, alpha = alpha, init = nmf.init, n.threads = n.cores, loss = loss)
  }

  colnames(nmf.res$W) <- rownames(nmf.res$H) <- sapply(1:ncol(nmf.res$W), function(i) paste("factor", i, sep = "_"))
  return(nmf.res)
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
#' @param alpha Regularization parameter
#' @param loss Loss function to use
#' @param n.cores Number of cores to use
#'
#' @return H matrix (scores) for newdata
#'
#' @import NNLM
#' @export
#'
ProjectSamples <- function(newdata, W, alpha = rep(0,3), loss = "mse", n.cores = 1) {
  lm.out <- NNLM::nnlm(W, newdata, alpha = alpha, loss = loss, n.threads = n.cores)
  return(lm.out$coefficients)
}

