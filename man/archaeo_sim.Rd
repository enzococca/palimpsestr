\name{archaeo_sim}
\alias{archaeo_sim}
\title{Simulate an archaeological palimpsest dataset}
\usage{archaeo_sim(n = 150, k = 3, seed = NULL, mixing = 0.08)}
\arguments{
  \item{n}{Number of observations.}
  \item{k}{Number of latent phases.}
  \item{seed}{Optional random seed.}
  \item{mixing}{Proportion of observations perturbed as mixed finds.}
}
\value{A data.frame with simulated archaeological observations.}
\description{Generate an exploratory synthetic dataset for testing the package workflow.}
