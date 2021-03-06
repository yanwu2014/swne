% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/similarity.R
\name{CalcKNN}
\alias{CalcKNN}
\title{kNN Graph Construction. Adapted from Seurat.}
\usage{
CalcKNN(data.use, k = 10, k.scale = 10)
}
\arguments{
\item{data.use}{Features x samples matrix to use to build the SNN}

\item{k}{Defines k for the k-nearest neighbor algorithm}

\item{k.scale}{Granularity option for k.param}
}
\value{
Returns kNN matrix in sparse matrix format
}
\description{
kNN Graph Construction. Adapted from Seurat.
}
