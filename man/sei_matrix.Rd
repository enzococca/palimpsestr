\name{sei_matrix}
\alias{sei_matrix}
\title{Compute the Stratigraphic Entanglement Index matrix}
\usage{sei_matrix(data, coords = c("x", "y", "z"), chrono = c("date_min", "date_max"), class_col = "class", weights = c(ws = 1, wz = 1, wt = 1, wc = 1), eps = 1e-09, z_floor = 0.25)}
\arguments{
  \item{data}{Input data.frame.}
  \item{coords}{Character vector of coordinate column names.}
  \item{chrono}{Character vector with minimum and maximum dating columns.}
  \item{class_col}{Class column name.}
  \item{weights}{Named numeric vector with SEI component weights.}
  \item{eps}{Small value to avoid division by zero.}
  \item{z_floor}{Minimum vertical denominator.}
}
\value{A symmetric numeric matrix.}
\description{Compute pairwise SEI values among observations.}
