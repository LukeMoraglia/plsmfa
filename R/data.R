#' Simulated data.
#'
#' This data set was created to demonstrate some of the properties of PLSMFA.
#'
#' @format A list of 6 items:
#' \describe{
#'   \item{\code{X}}{A matrix of 100 observations and 9 variables, including two submatrices}
#'   \item{\code{Y}}{A matrix of 100 observations and 6 variables, including two submatrices}
#'   \item{\code{X_design}}{A vector indicating which submatrix of X each column belongs too}
#'   \item{\code{Y_design}}{A vector indicating which submatrix of Y each column belongs too}
#'   \item{\code{X_var_color}}{A list containing color vectors for the variables, \code{$oc},
#'    and for the submatrices, \code{$gc}, of X}
#'   \item{\code{Y_var_color}}{A list containing color vectors for the variables, \code{$oc},
#'    and for the submatrices, \code{$gc}, of Y}
#'   \item{\code{obs_design}}{A length 100 vector indicating the (fictional) group of each observation.}
#'   \item{\code{obs_color}}{Colors for the observations based on group}
#'   \item{\code{group_color}}{Colors for the 3 groups of observations}
#' }
"sim_data"
