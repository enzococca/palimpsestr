\name{ese}
\alias{ese}
\title{Compute Excavation Stratigraphic Energy}
\usage{ese(data, coords = c("x", "y", "z"), chrono = c("date_min", "date_max"), class_col = "class", beta = c(1, 1, 1, 1), neighbourhood = NULL)}
\arguments{
  \item{data}{Input data.frame.}
  \item{coords}{Character vector of coordinate column names.}
  \item{chrono}{Character vector with minimum and maximum dating columns.}
  \item{class_col}{Class column name.}
  \item{beta}{Weights for spatial, vertical, temporal and class mismatch.}
  \item{neighbourhood}{Optional maximum XY distance.}
}
\value{A numeric vector with local energy values.}
\description{Measure local depositional tension or incompatibility among nearby observations.}
