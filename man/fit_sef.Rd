\name{fit_sef}
\alias{fit_sef}
\title{Fit the Stratigraphic Entanglement Field model}
\usage{fit_sef(data, coords = c("x", "y", "z"), chrono = c("date_min", "date_max"), class = "class", tafonomy = NULL, context = NULL, harris = NULL, k = 3, weights = c(ws = 1, wz = 1, wt = 1, wc = 1), seed = 1, em_iter = 10)}
\arguments{
  \item{data}{Input data.frame.}
  \item{coords}{Character vector of coordinate column names.}
  \item{chrono}{Character vector with minimum and maximum dating columns.}
  \item{class}{Class/material column name.}
  \item{tafonomy}{Optional tafonomy score column.}
  \item{context}{Optional context/SU column.}
  \item{harris}{Optional n x n matrix of pairwise penalties.}
  \item{k}{Number of phases.}
  \item{weights}{Named numeric vector used by \code{sei_matrix()}.}
  \item{seed}{Optional random seed.}
  \item{em_iter}{Number of soft reallocation iterations.}
}
\value{An object of class \code{sef_fit}.}
\description{Fit an exploratory SEF model with soft phase assignment and diagnostic outputs.}
