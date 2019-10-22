#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <math.h>
#include <vector>
#include <stdexcept>
#include <queue>
#include <cmath>
#include <unordered_map>
// #include <fstream>
// #include <string>
// #include <iomanip>
// include "data_manipulation.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// Adapted from pagoda2: https://github.com/hms-dbmi/pagoda2
// calculates factor-stratified sums for each column
// rowSel is an integer factor;
// note that the 0-th column will return sums for any NA values; 0 or negative values will be omitted
// [[Rcpp::export]]
arma::mat colSumByFac(SEXP sY,  SEXP rowSel) {
// need to do this as SEXP, modify the slots on the fly
  S4 mat(sY);
  const arma::uvec i(( unsigned int *)INTEGER(mat.slot("i")), LENGTH(mat.slot("i")), false, true);
  const arma::ivec dims(INTEGER(mat.slot("Dim")), LENGTH(mat.slot("Dim")), false, true);
  const arma::ivec p(INTEGER(mat.slot("p")), LENGTH(mat.slot("p")), false, true);
  arma::vec Y(REAL(mat.slot("x")), LENGTH(mat.slot("x")), false, true);

  const arma::ivec rs = arma::ivec(INTEGER(rowSel), LENGTH(rowSel), false, true);

  int ncols = p.size() - 1;
  int nlevels = 0;
  for(int j = 0; j < rs.size(); j++) {
    if(rs[j] != NA_INTEGER) {
      if(rs[j] > nlevels) { nlevels = rs[j]; }
    }
  }
  if(nlevels == 0) { stop("colSumByFac(): supplied factor doesn't have any levels!"); }
  arma::mat sumM(nlevels + 1, ncols, arma::fill::zeros);

  // for each gene
  // pragma omp parallel for shared(sumM)
  for(int g = 0; g < ncols; g++) {
    int p0 = p[g]; int p1 = p[g + 1];
    if(p1 - p0 < 1) { continue; }
    for(int j = p0; j < p1; j++) {
      int row = i[j];
      int f = rs[row];
      if(f == NA_INTEGER) {
	      sumM(0, g) += Y[j];
      } else if(f > 0) {
      	sumM(f, g) += Y[j];
      }
    }
  }
  return sumM;
}

// Adapted from pagoda2: https://github.com/hms-dbmi/pagoda2
// calculate column mean and variance, optionally taking a subset of rows to operate on
// [[Rcpp::export]]
Rcpp::DataFrame colMeanVarS(SEXP sY,  SEXP rowSel) {
// need to do this as SEXP, modify the slots on the fly
  S4 mat(sY);
  const arma::uvec i(( unsigned int *)INTEGER(mat.slot("i")),LENGTH(mat.slot("i")),false,true);
  const arma::ivec dims(INTEGER(mat.slot("Dim")),LENGTH(mat.slot("Dim")),false,true);
  const arma::ivec p(INTEGER(mat.slot("p")),LENGTH(mat.slot("p")),false,true);
  arma::vec Y(REAL(mat.slot("x")),LENGTH(mat.slot("x")),false,true);

  bool rowSelSpecified=!Rf_isNull(rowSel);
  const arma::ivec rs=(rowSelSpecified) ? arma::ivec(INTEGER(rowSel),LENGTH(rowSel),false,true) : arma::ivec();

  int ncols=p.size()-1;
  int nrows=dims[0];
  if(rowSelSpecified) {
    nrows=0;
    for(int j=0;j<rs.size();j++) { if(rs[j]) { nrows++; } }
  }
  arma::vec meanV(ncols,arma::fill::zeros); arma::vec varV(ncols,arma::fill::zeros); arma::vec nobsV(ncols,arma::fill::zeros);
  // for each gene
#pragma omp parallel for shared(meanV,varV,nobsV)
  for(int g=0;g<ncols;g++) {
    int p0=p[g]; int p1=p[g+1];
    if(p1-p0 <1) { continue; }
    arma::colvec ly;
    if(rowSelSpecified) {
      // select valid rows
      int nvalid=0;
      ly=arma::vec(p1-p0);
      for(int j=p0;j<p1;j++) {
	if(rs[i[j]]) {
	  ly[nvalid]=Y[j]; nvalid++;
	}
      }
      nobsV[g]=nvalid;
      ly=ly.head(nvalid);
    } else {
      nobsV[g]=p1-p0;
      ly=Y.subvec(p0,p1-1);
    }

    double m=sum(ly)/nrows;
    meanV[g]=m;
    ly-=m; ly%=ly;
    varV[g]=(sum(ly)+(m*m*(nrows-ly.size())))/nrows;
  }
  return Rcpp::DataFrame::create(Named("m")=meanV,Named("v")=varV,Named("nobs",nobsV));
}

// Adapted from pagoda2: https://github.com/hms-dbmi/pagoda2
// Winsorize top N values in each column of a sparse matrix
// [[Rcpp::export]]
int inplaceWinsorizeSparseCols(SEXP sY, const int n) {
  // need to do this as SEXP, modify the slots on the fly
  S4 mat(sY);
  const arma::uvec i(( unsigned int *)INTEGER(mat.slot("i")),LENGTH(mat.slot("i")),false,true);
  const arma::ivec dims(INTEGER(mat.slot("Dim")),LENGTH(mat.slot("Dim")),false,true);
  const arma::ivec p(INTEGER(mat.slot("p")),LENGTH(mat.slot("p")),false,true);
  arma::vec Y(REAL(mat.slot("x")),LENGTH(mat.slot("x")),false,true);
  int ncols=p.size()-1;
  // for each column

  arma::vec tv(ncols);

  for(int g=0;g<ncols;g++) {
    std::priority_queue<std::pair<double,int> , std::vector<std::pair<double,int> >, std::greater<std::pair<double,int> > > q;
    int p1=p[g+1]; int p0=p[g]; int ncells=p1-p0;
    if(ncells<=(n+1)) { continue; } // not enough observations.. could 0 them all out, alternatively, but presumably this was done beforehand
    // find indices of top n+1 values
    // insert first n+1 values
    for(int j=p0;j<p0+n+1;j++) {
      q.push(std::pair<double,int>(Y[j],j));
    }
    for(int j=p0+n+1;j<p1;j++) {
      double v=Y[j];
      if(v>q.top().first) {
        q.pop();
        q.push(std::pair<double,int>(v,j));
      }
    }
    // set all values to the smallest (top) one
    double v=q.top().first;
    q.pop();
    while(!q.empty()) {
      Y[q.top().second]=v;
      q.pop();
    }
  }
  return(1);
}


typedef Eigen::Triplet<double> T;
// Adapted from Seurat
//[[Rcpp::export]]
Eigen::SparseMatrix<double> ComputeSNN(Eigen::MatrixXd nn_ranked, double prune) {
  std::vector<T> tripletList;
  int k = nn_ranked.cols();
  tripletList.reserve(nn_ranked.rows() * nn_ranked.cols());
  for(int j=0; j<nn_ranked.cols(); ++j){
    for(int i=0; i<nn_ranked.rows(); ++i) {
      tripletList.push_back(T(i, nn_ranked(i, j) - 1, 1));
    }
  }
  Eigen::SparseMatrix<double> SNN(nn_ranked.rows(), nn_ranked.rows());
  SNN.setFromTriplets(tripletList.begin(), tripletList.end());
  SNN = SNN * (SNN.transpose());
  for (int i=0; i < SNN.outerSize(); ++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it){
      it.valueRef() = it.value()/(k + (k - it.value()));
      if(it.value() < prune){
        it.valueRef() = 0;
      }
    }
  }
  SNN.prune(0.0); // actually remove pruned values
  return SNN;
}


// Compute statistical significance of number of shared edges between two clusters
// Adapted from monocle3 (https://github.com/cole-trapnell-lab/monocle3/blob/master/src/clustering.cpp)
NumericMatrix pnorm_over_mat_cpp(NumericMatrix num_links_ij, NumericMatrix var_null_num_links) {
  int n = num_links_ij.nrow();
  NumericMatrix tmp(n, n);

  for (int i = 0; i < n; i ++) {
    for (int j = 0; j < n; j ++) {
      // tmp(i, j) = Rcpp::pnorm( num_links_ij(i, j), 0.0, sqrt(var_null_num_links(i, j)), bool lower = false, bool log = false );
      tmp(i, j) = R::pnorm(num_links_ij(i, j), 0.0, sqrt(var_null_num_links(i, j)), 0, 0);
    }
  }
  return tmp;
}


// Compute statistical significance of number of shared edges between two clusters
// Adapted from monocle3 (https://github.com/cole-trapnell-lab/monocle3/blob/master/src/clustering.cpp)
// [[Rcpp::export]]
NumericMatrix pnorm_over_mat(SEXP R_num_links_ij, SEXP R_var_null_num_links) {
  NumericMatrix num_links_ij(R_num_links_ij);
  NumericMatrix var_null_num_links(R_var_null_num_links);

  return pnorm_over_mat_cpp(num_links_ij, var_null_num_links);
}
