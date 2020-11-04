#' Compute partial latent variables for a subtable
#'
#' Whereas latent variables are computed as linear combinations of all the variables in a data
#' table, partial latent variables are computed from only the variables of a specific subtable.
#'
#'
#' @param data_ori Should be the MFA normalized data, so either \code{mfa_norm_X} or \code{mfa_norm_Y}
#' @param loading_mat The corresponding matrix of loadings, either \code{u} or \code{v}
#' @param from_index The starting index of the subtable
#' @param to_index The ending index of the subtable
#' @param n_subtable The total number of subtables in \code{data_ori}
#'
#' @return The partial latent variables for the subtable
#' @author Luke Moraglia
partial_latent_variable <- function(data_ori,
                                    loading_mat,
                                    from_index,
                                    to_index,
                                    n_subtable) {
   data_ori <- as.matrix(data_ori)
   loading_mat <- as.matrix(loading_mat)

   plv <- n_subtable *
      (data_ori[,from_index:to_index] %*%
          loading_mat[from_index:to_index,])

   return(plv)
}


#' Compute the group means of the partial latent variables, based on a design vector
#'
#' @param partial_lv A block of partial latent variables. Should be N x dims x subtables
#' @param design A design vector for the N observations
#'
#' @return A block of plv means. Groups x dims x subtables
#' @author Luke Moraglia
get_partial_lv_group_means <- function(partial_lv,
                                       design){
   partial_lv_means <- purrr::map(1:dim(partial_lv)[3],
                                  ~as.matrix(PTCA4CATA::getMeans(partial_lv[,,.], design))) %>%
      simplify2array()
   dimnames(partial_lv_means)[[3]] <- dimnames(partial_lv)[[3]]
   return(partial_lv_means)
}
