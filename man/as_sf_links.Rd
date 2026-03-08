\name{as_sf_links}
\alias{as_sf_links}
\title{Export high-SEI links as an sf LINESTRING layer}
\usage{
as_sf_links(object, quantile_threshold = 0.9, crs = NA)
}
\arguments{
\item{object}{A \code{sef_fit} object.}
\item{quantile_threshold}{Quantile used to retain only the strongest links.}
\item{crs}{Optional CRS for the output geometry.}
}
\value{
An \code{sf} object with one line per retained pair.
}
\description{
Export high-SEI links as an sf LINESTRING layer.
}
