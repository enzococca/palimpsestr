\name{as_sf_phase}
\alias{as_sf_phase}
\title{Convert fitted results to sf}
\usage{as_sf_phase(object, crs = NA)}
\arguments{\item{object}{A \code{sef_fit} object.}
\item{crs}{Optional CRS passed to \code{sf::st_as_sf()}.}}
\value{An \code{sf} object.}
\description{Converts fitted phase diagnostics to an \code{sf} point layer when the \code{sf} package is installed.}
