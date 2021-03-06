% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plsmfa_graphs.R
\name{correlation_circle}
\alias{correlation_circle}
\title{Create a correlation circle plot for either X or Y}
\usage{
correlation_circle(
  res_plsmfa,
  axis1 = 1,
  axis2 = 2,
  XorY = "X",
  col_vars = "blueviolet"
)
}
\arguments{
\item{res_plsmfa}{output of \code{\link{plsmfa}}}

\item{axis1}{(default 1) the latent variable to plot on the horizontal axis}

\item{axis2}{(default 2) the latent variable to plot on the vertical axis}

\item{XorY}{(default "X") Either "X" or "Y", showing which matrix the LVs should come from}

\item{col_vars}{(default "blueviolet") Either one color or a vector of colors
with a length of N variables}
}
\value{
A list containing:
 \item{\code{plot}}{The final plot from the function}
 \item{\code{base_plot}}{Output from \code{\link[PTCA4CATA]{createFactorMap}}}
 \item{\code{arrows}}{The output of \code{\link[PTCA4CATA]{addArrows}}}
 \item{\code{circle}}{The output of \code{\link[PTCA4CATA]{addCircleOfCor}}}
}
\description{
Create a correlation circle plot for either X or Y
}
\author{
Luke Moraglia
}
