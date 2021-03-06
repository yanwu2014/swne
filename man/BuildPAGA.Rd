% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paga.R
\name{BuildPAGA}
\alias{BuildPAGA}
\title{Make PAGA graph. Adapted from Monocle3}
\usage{
BuildPAGA(knn, clusters = NULL, qval.cutoff = 0.05)
}
\arguments{
\item{knn}{kNN in the form of a sparse matrix}

\item{clusters}{Cell clusters}

\item{qval.cutoff}{Q-value cutoff for cluster to cluster interactions}
}
\value{
List representing a PAGA graph
}
\description{
Make PAGA graph. Adapted from Monocle3
}
