#' Perform partial least squares multiple factor analysis (PLSMFA)
#'
#' PLSMFA is a generalization of PLSC when the two data tables have subtables of variables.
#' The matrices \code{data1} and \code{data2} must have columns that are ordered by groups.
#' The column group memberships are stored in \code{column_design1} and \code{column_design2}.
#' If your data's columns are not ordered by their groups, but you have the out of order column
#' design vector, you can use \code{\link{arrange_data_by_design}} to order the columns.
#'
#' PLSMFA implements PLSC with an additional preprocessing step that normalizes subtables within
#' the two data tables to be analyzed. The normalization is borrowed from multiple factor analysis (MFA),
#' which normalizes subtables by dividing them by their first singular value.
#'
#' PLSMFA begins by centering
#' and scaling the variables as given in \code{center1}, \code{center2}, \code{scale1}, and \code{scale2}.
#' Then each subtable is decomposed by the SVD and divided by its first singular value (\code{\link{mfa_norm}}).
#' The cross-product matrix of the normalized data tables is computed and decomposed by the SVD.
#'
#' The results are very similar to PLSC. However, \code{plsmfa} also computes the contribution of
#' each variable and each subtable. Additionally, each subtable has a corresponding matrix of partial latent
#' variables, which describe the observations as seen by only this subtable.
#'
#'
#' @param data1 The first data matrix, also called "X", with \eqn{N} rows and \eqn{I} columns
#'              organized into \eqn{K} subtables.
#' @param data2 The second data matrix, also called "Y", with \eqn{N} rows and \eqn{J} columns
#'              organized into \eqn{M} subtables.
#' @param column_design1 A vector of length \eqn{I} that gives the subtable membership of the X variables
#' @param column_design2 A vector of length \eqn{J} that gives the subtable membership of the Y variables
#' @param center1 (default TRUE) Whether to center \code{data1}
#' @param center2 (default TRUE) Whether to center \code{data2}
#' @param scale1 (default "SS1") Whether to scale \code{data1}.
#'  Any acceptable input to \code{\link[ExPosition]{expo.scale}}
#' @param scale2 (default "SS1") Whether to scale \code{data2}
#'  Any acceptable input to \code{\link[ExPosition]{expo.scale}}
#'
#' @return A list containing:
#'   \item{\code{pls}}{results from the underlying PLSC analysis. This is output from \code{\link[GSVD]{gplssvd}}.
#'   It is a list containing:
#'    \describe{
#'     \item{\code{d}}{The singular values}
#'     \item{\code{u}}{The left singular vectors. For PLSMFA, these are the loadings of X.}
#'     \item{\code{v}}{The right singular vectors. For PLSMFA, these are the loadings of Y.}
#'     \item{\code{d_full}}{Should be the same as \code{d}}
#'     \item{\code{l_full}}{Should be the same as \code{l}}
#'     \item{\code{l}}{The eigenvalues}
#'     \item{\code{lx}}{The latent variables from X}
#'     \item{\code{ly}}{The latent variables from Y}
#'     \item{\code{p}}{The left generalized singular vectors. Should be the same as \code{u}}
#'     \item{\code{fi}}{The X variable "factor scores".
#'                      This is the left singular vectors multiplied by their singular values.}
#'     \item{\code{q}}{The right generalized singular vectors. Should be the same as \code{v}}
#'     \item{\code{fj}}{The Y variable "factor scores".
#'                      This is the right singular vectors multiplied by their singular values.}
#'     }
#'   }
#'   \item{\code{normed_X}}{X after MFA normalization of the subtables}
#'   \item{\code{normed_Y}}{Y after MFA normalization of the subtables}
#'   \item{\code{first_singval_X}}{A vector containing the first singular value of each subtable in X}
#'   \item{\code{first_singval_Y}}{A vector containing the first singular value of each subtable in Y}
#'   \item{\code{partial_lv_X}}{A 3D array, of size \eqn{N \times L \text{ (number of dimensions) } \times K}
#'    (subtables of X), containing the partial latent variables for the X subtables}
#'   \item{\code{partial_lv_Y}}{A 3D array, of size \eqn{N \times L \text{ (number of dimensions) } \times M}
#'    (subtables of Y), containing the partial latent variables for the Y subtables}
#'   \item{\code{design}}{A list containing information about the column designs including:
#'                 \describe{
#'                 \item{\code{column_design_X}}{the same as \code{column_design1} from input}
#'                 \item{\code{column_design_Y}}{the same as \code{column_design2} from input}
#'                 \item{\code{from_indices_X}}{a vector containing the column index that starts each subtable of X}
#'                 \item{\code{to_indices_X}}{a vector containing the column index that ends each subtable of X}
#'                 \item{\code{from_indices_Y}}{a vector containing the column index that starts each subtable of Y}
#'                 \item{\code{to_indices_Y}}{a vector containing the column index that ends each subtable of Y}
#'                 }
#'                 }
#'    \item{\code{ctr}}{Contributions of variables and subtables. Output from \code{\link{contributions}}}
#'
#' @author Luke Moraglia
#' @export
#'
#' @importFrom magrittr %>%
plsmfa <- function(data1, data2,
                     column_design1, column_design2,
                     center1 = TRUE, center2 = TRUE,
                     scale1 = "SS1", scale2 = "SS1"){

   # Center and scale
   ZX <- ExPosition::expo.scale(data1, center1, scale1)
   ZY <- ExPosition::expo.scale(data2, center2, scale2)

   # column designs become group by variable indicator matrices
   column_design_X <- t(ExPosition::makeNominalData(as.matrix(column_design1)))
   column_design_Y <- t(ExPosition::makeNominalData(as.matrix(column_design2)))

   design_row_sum_X <- rowSums(column_design_X)
   design_row_sum_Y <- rowSums(column_design_Y)

   # to and from indices show the start and end index of each group of columns
   to_indices_X <- purrr::accumulate(design_row_sum_X, `+`)
   to_indices_Y <- purrr::accumulate(design_row_sum_Y, `+`)

   from_indices_X <- dplyr::lag(to_indices_X, n = 1, default = 0) + 1
   from_indices_Y <- dplyr::lag(to_indices_Y, n = 1, default = 0) + 1

   names(from_indices_X) <- names(to_indices_X)
   names(from_indices_Y) <- names(to_indices_Y)

   n_subtable_X <- length(from_indices_X)
   n_subtable_Y <- length(from_indices_Y)

   # mfa_norm_res is a list containing the normalized matrices and the first singular values
   mfa_norm_res_X <- purrr::map2(from_indices_X, to_indices_X, ~mfa_norm(ZX, .x, .y))
   mfa_norm_res_Y <- purrr::map2(from_indices_Y, to_indices_Y, ~mfa_norm(ZY, .x, .y))

   # grabs the normalized data
   mfa_norm_X <- purrr::map(mfa_norm_res_X, 1) %>% purrr::reduce(cbind)
   mfa_norm_Y <- purrr::map(mfa_norm_res_Y, 1) %>% purrr::reduce(cbind)

   # grabs the first singular values
   first_singval_X <- purrr::map_dbl(mfa_norm_res_X, 2)
   first_singval_Y <- purrr::map_dbl(mfa_norm_res_Y, 2)

   colnames(mfa_norm_X) <- colnames(ZX)
   colnames(mfa_norm_Y) <- colnames(ZY)

   # the underlying PLSC call. first_pos_gplssvd makes the results consistent across machines
   res_pls <- first_pos_gplssvd(GSVD::gplssvd(mfa_norm_X, mfa_norm_Y), mfa_norm_X, mfa_norm_Y)

   # Partial latent variables
   partial_lv_X <- purrr::map2(from_indices_X, to_indices_X,
                               ~partial_latent_variable(mfa_norm_X,
                                                        res_pls$u,
                                                        .x, .y,
                                                        n_subtable_X)) %>%
      simplify2array()

   partial_lv_Y <- purrr::map2(from_indices_Y, to_indices_Y,
                               ~partial_latent_variable(mfa_norm_Y,
                                                        res_pls$v,
                                                        .x, .y,
                                                        n_subtable_Y)) %>%
      simplify2array()

   ctr <- contributions(res_pls$u, res_pls$v, column_design1, column_design2)

   design_list <- list(column_design_X = column_design1,
                       column_design_Y = column_design2,
                       from_indices_X = from_indices_X,
                       to_indices_X = to_indices_X,
                       from_indices_Y = from_indices_Y,
                       to_indices_Y = to_indices_Y)

   res_plsmfa <- list(pls = res_pls,
                       normed_X = mfa_norm_X,
                       normed_Y = mfa_norm_Y,
                       first_singval_X = first_singval_X,
                       first_singval_Y = first_singval_Y,
                       partial_lv_X = partial_lv_X,
                       partial_lv_Y = partial_lv_Y,
                       design = design_list,
                      ctr = ctr
                      )

   invisible(res_plsmfa)
}

#' Normalize a matrix by its first singular value
#'
#' Divides a matrix by its first singular value, the norm used in MFA and PLSMFA.
#'
#' @param data The matrix to be normalized
#' @param from (default 1) If subsetting the columns of \code{data}, the starting index
#' @param to (default \code{NCOL(data)}) If subsetting the columns of \code{data}, the ending index
#'
#' @return A list containing
#'   \enumerate{
#'   \item{The normalized data}
#'   \item{The first singular value used to normalize}
#'   }
#' @author Luke Moraglia
#' @export
#'
mfa_norm <- function(data, from = 1, to = NCOL(data)){
   sing_1 <- svd(data[,from:to], nu = 0, nv = 0)$d[1]
   list(data[,from:to]/sing_1,
        sing_1)
}


#' Make the first element of each singular vector positive for results of \code{GSVD::gplssvd}
#'
#' This is a helper function that prevents problems related to axis flipping by setting the first
#' element of each singular vector in \code{u} as positive.
#'
#' Each column of \code{u} is set so that the first element is positive. The rest of the results
#' are adapted accordingly.
#'
#' @param res_gplssvd Output from \code{\link[GSVD]{gplssvd}}
#' @param X Input to \code{gplssvd}
#' @param Y Input to \code{gplssvd}
#'
#' @return Results from gplssvd changed accordingly.
#' @author Luke Moraglia
first_pos_gplssvd <- function(res_gplssvd, X, Y){
   res_firstpos <- data4PCCAR::firstpos(res_gplssvd$u, res_gplssvd$v)
   res_gplssvd$u <- res_firstpos$P
   res_gplssvd$v <- res_firstpos$Q
   res_gplssvd$lx <- X %*% res_gplssvd$u
   res_gplssvd$ly <- Y %*% res_gplssvd$v
   res_firstpos <- data4PCCAR::firstpos(res_gplssvd$p, res_gplssvd$q)
   res_gplssvd$p <- res_firstpos$P
   res_gplssvd$q <- res_firstpos$Q
   res_firstpos <- data4PCCAR::firstpos(res_gplssvd$fi, res_gplssvd$fj)
   res_gplssvd$fi <- res_firstpos$P
   res_gplssvd$fj <- res_firstpos$Q
   return(res_gplssvd)
}
