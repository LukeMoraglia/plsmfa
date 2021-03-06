% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plsmfa_graphs.R
\name{partial_lv_group_means_map}
\alias{partial_lv_group_means_map}
\title{Plot partial latent variable group means}
\usage{
partial_lv_group_means_map(
  res_plsmfa,
  design,
  lv_plot,
  axis1 = 1,
  axis2 = 2,
  XorY = "X",
  col_items = NULL,
  col_blocks = NULL,
  use_defaults = TRUE,
  ...
)
}
\arguments{
\item{res_plsmfa}{output of \code{\link{plsmfa}}}

\item{design}{A design vector for the rows}

\item{lv_plot}{the output of \code{\link{latent_variable_map}}}

\item{axis1}{(default 1) the latent variable to plot on the horizontal axis}

\item{axis2}{(default 2) the latent variable to plot on the vertical axis}

\item{XorY}{(default "X") Either "X" or "Y", showing which matrix the LVs should come from}

\item{col_items}{(default NULL)  Either one color or a vector of colors with a length of N groups}

\item{col_blocks}{(default NULL)  Either one color or a vector of colors with a length of N subtables}

\item{use_defaults}{(default TRUE) set to FALSE to use other \code{...} args}

\item{...}{args passed to \code{\link[PTCA4CATA]{createPartialFactorScoresMap}}}
}
\value{
A list containing:
  \item{\code{plot}}{The final plot from the function}
  \item{\code{part_lv_means_plot}}{Output of \code{\link[PTCA4CATA]{createPartialFactorScoresMap}}}
  \item{\code{lv_plot}}{The \code{lv_plot} from input}
}
\description{
Plot partial latent variable group means
}
\author{
Luke Moraglia
}
