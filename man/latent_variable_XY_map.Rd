% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plsmfa_graphs.R
\name{latent_variable_XY_map}
\alias{latent_variable_XY_map}
\title{Latent variable map for a single dimension, X vs Y}
\usage{
latent_variable_XY_map(
  res_plsmfa,
  lv_num = 1,
  design = NULL,
  col_obs = "blueviolet",
  col_group = "blueviolet",
  n_iter = 100,
  xlab = paste0("Latent Variable X", lv_num),
  ylab = paste0("Latent Variable Y", lv_num),
  ...
)
}
\arguments{
\item{res_plsmfa}{output of \code{\link{plsmfa}}}

\item{lv_num}{(default 1) which dimension's latent variables to plot}

\item{design}{(default NULL) A design vector for the rows. If included, group means are plotted}

\item{col_obs}{(default "blueviolet") Either one color or a vector of colors with a length of N observations}

\item{col_group}{(default "blueviolet") Either one color or a vector of colors with a length of N groups}

\item{n_iter}{(default 100) Number of bootstrap iterations for computing confidence intervals around the group means}

\item{xlab}{(default \code{paste0("Latent Variable X", lv_num)}) Horizontal axis label}

\item{ylab}{(default \code{paste0("Latent Variable Y", lv_num)}) Vertical axis label}

\item{...}{args passed to \code{\link[PTCA4CATA]{createFactorMap}} for the lv_fact_map}
}
\value{
A list containing:
  \item{\code{plot}}{The final plot from the function}
  \item{\code{lv_fact_map}}{Output from \code{\link[PTCA4CATA]{createFactorMap}} for the observations}
  \item{\code{group_means}}{The group means on the latent variables}
  \item{\code{lv_group_means_fact_map}}{Output from \code{\link[PTCA4CATA]{createFactorMap}} for the group means}
  \item{\code{CI_ellipses}}{The confidence intervals of the group means}
}
\description{
Latent variable map for a single dimension, X vs Y
}
\author{
Luke Moraglia
}
