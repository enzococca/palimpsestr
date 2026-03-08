\name{compare_k}
\alias{compare_k}
\title{Compare multiple candidate phase counts}
\usage{compare_k(data, k_values = 2:6, ...)}
\arguments{\item{data}{Input data.frame.}
\item{k_values}{Integer vector of candidate phase counts.}
\item{...}{Additional arguments passed to \code{fit_sef()}.}}
\value{A data.frame with exploratory fit diagnostics.}
\description{Fits the prototype model across several values of \code{k} and returns a compact comparison table.}
