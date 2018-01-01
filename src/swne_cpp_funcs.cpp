// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <vector>
#include <stdexcept>
#include <queue>
using namespace Rcpp;


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
