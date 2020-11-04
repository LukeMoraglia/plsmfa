#' Permutation tests for the eigenvalues of PLSMFA
#'
#' The permutation test randomizes each column of the (normalized) data tables and performs PLSMFA
#' on the permuted data tables. This randomization simulates the pattern expected under the null hypothesis of no relationships
#' between the variables. This process repeats \code{n_iter} times to create a null distribution of the eigenvalues, from
#' which probabilities for the eigenvalues, total inertia, and percentage of inertia are computed.
#'
#' Some additional details: PLSMFA starts off by normalizing each subtable, which changes the inertia of each subtable.
#' If the original data were permuted, the first singular value of each subtable would change, and in the case of the
#' unidimensional subtables, the first singular value would decrease, sometimes by a lot. When the first singular value
#' decreases, the subtable is not reduced by as much, and so the normalized permuted subtable would have more inertia than the
#' normalized original subtable. This translates into more inertia for the cross-product matrix that is decomposed by the
#' SVD, and so it would make comparisons between the eigenvalues of the original result and the permuted results difficult.
#'
#' The solution implemented here is to instead permute the normalized data tables, so that there are no relationships between
#' the variables of the normalized data tables. However, this means that the permutation test acts as if the first singular
#' value of each subtable is a fixed value, which may not necessarily be the case. In order to add variability to the
#' first singular value, but still keep the inertia comparable between iterations, there is the option to bootstrap the
#' first singular values. This option (\code{bootstrap_first_singval}) is similar to the process described in
#' \code{\link{boot_plsmfa}} where samples are drawn with replacement to simulate fluctuations in the data due to sampling.
#' While this option may increase the robustness of the permutation test, it also adds considerable computational time
#' which could be apparent when the data have many variables and/or subtables.
#'
#' @inheritParams plsmfa
#' @param n_iter (default 1000) The number of permutations to compute
#' @param compact (default FALSE) If TRUE, returns only the p values. If FALSE, includes the results of each iteration.
#' @param bootstrap_first_singval (default FALSE) If TRUE, uses a bootstrap sample
#'   to compute the first singular values of the subtables. Setting to TRUE adds
#'   considerable computation time.
#'
#' @return A list containing:
#'    \item{\code{p_eig}}{the p value of each eigenvalue}
#'    \item{\code{p_inertia}}{the p value for the total inertia}
#'    \item{\code{p_percent_inertia}}{the p values for the percentage of inertia explained by each dimension}
#'    and if \code{compact = FALSE}:
#'    \item{\code{fixed_eig}}{the observed eigenvalues from PLSMFA}
#'    \item{\code{fixed_inertia}}{the observed inertia from PLSMFA}
#'    \item{\code{fixed_percent_inertia}}{the obeserved percentage of inertia explained by each dimension from PLSMFA}
#'    \item{\code{perm_eig}}{the eigenvalues from each permutation}
#'    \item{\code{perm_inertia}}{the inertia from each permutation}
#'    \item{\code{perm_percent_inertia}}{the percentage of inertia explained by the dimensions of each permutation}
#' @author Luke Moraglia
#' @export
perm_plsmfa <- function(data1, data2,
                        column_design1, column_design2,
                        center1 = TRUE, center2 = TRUE,
                        scale1 = "SS1", scale2 = "SS1",
                        n_iter = 1000,
                        compact = FALSE,
                        bootstrap_first_singval = FALSE){
   ZX <- ExPosition::expo.scale(data1, center1, scale1)
   ZY <- ExPosition::expo.scale(data2, center2, scale2)

   column_design_X <- t(ExPosition::makeNominalData(column_design1))
   column_design_Y <- t(ExPosition::makeNominalData(column_design2))

   design_row_sum_X <- rowSums(column_design_X)
   design_row_sum_Y <- rowSums(column_design_Y)

   to_indices_X <- purrr::accumulate(design_row_sum_X, `+`)
   to_indices_Y <- purrr::accumulate(design_row_sum_Y, `+`)

   from_indices_X <- dplyr::lag(to_indices_X, n = 1, default = 0) + 1
   from_indices_Y <- dplyr::lag(to_indices_Y, n = 1, default = 0) + 1

   names(from_indices_X) <- names(to_indices_X)
   names(from_indices_Y) <- names(to_indices_Y)

   n_subtable_X <- length(from_indices_X)
   n_subtable_Y <- length(from_indices_Y)

   mfa_norm_res_X <- purrr::map2(from_indices_X, to_indices_X, ~mfa_norm(ZX, .x, .y))
   mfa_norm_res_Y <- purrr::map2(from_indices_Y, to_indices_Y, ~mfa_norm(ZY, .x, .y))

   mfa_norm_X <- purrr::map(mfa_norm_res_X, 1) %>% purrr::reduce(cbind)
   mfa_norm_Y <- purrr::map(mfa_norm_res_Y, 1) %>% purrr::reduce(cbind)

   colnames(mfa_norm_X) <- colnames(ZX)
   colnames(mfa_norm_Y) <- colnames(ZY)

   # The fixed results from PLSMFA
   fixed_singval <- GSVD::tolerance_svd(t(mfa_norm_X) %*% mfa_norm_Y, nu = 1, nv = 1)$d
   fixed_eig <- fixed_singval^2
   fixed_inertia <- sum(fixed_eig)
   fixed_percent_inertia <- fixed_eig / fixed_inertia
   n_singval <- length(fixed_singval)


   if(bootstrap_first_singval){
      perm_eig <- replicate(n_iter, .permute_bootstrap_first_singval(data1, data2, n_singval + 1,
                                                                  from_indices_X, to_indices_X,
                                                                  from_indices_Y, to_indices_Y,
                                                                  center1, center2,
                                                                  scale1, scale2))
   } else {
      perm_eig <- replicate(n_iter,
                            .permute_fixed_first_singval(mfa_norm_X, mfa_norm_Y, n_singval + 1))
   }


   perm_eig <- t(perm_eig[1:n_singval, ])
   perm_inertia <- rowSums(perm_eig)
   perm_percent_inertia <- perm_eig / perm_inertia

   p_eig <- rowSums(t(perm_eig) > fixed_eig) / n_iter
   p_eig[p_eig == 0] <- 1 / n_iter
   p_inertia <- sum(perm_inertia > fixed_inertia) / n_iter
   if(p_inertia == 0){
      p_inertia <- 1 / n_iter
   }
   p_percent_inertia <- rowSums(t(perm_percent_inertia) > fixed_percent_inertia) / n_iter
   p_percent_inertia[p_percent_inertia == 0] <- 1 / n_iter

   return_list <- list(p_eig = p_eig,
                       p_inertia = p_inertia,
                       p_percent_inertia = p_percent_inertia)
   if(!compact){
      return_list$fixed_eig = fixed_eig
      return_list$fixed_inertia = fixed_inertia
      return_list$fixed_percent_inertia = fixed_percent_inertia
      return_list$perm_eig = perm_eig
      return_list$perm_inertia = perm_inertia
      return_list$perm_percent_inertia = perm_percent_inertia
   }


   return(return_list)
}


.permute_fixed_first_singval <- function(X, Y, length){
   X_rand <- apply(X, 2, sample)
   Y_rand <- apply(Y, 2, sample)
   R_rand <- t(X_rand) %*% Y_rand

   eig_rand <- rep(0, length)
   # This is primarily to increase speed. eigen is faster than svd.
   if(dim(X)[2] < dim(Y)[2]){
      res_eig <- eigen(R_rand %*% t(R_rand), symmetric = TRUE, only.values = TRUE)$values
   } else {
      res_eig <- eigen(t(R_rand) %*% R_rand, symmetric = TRUE, only.values = TRUE)$values
   }

   eig_rand[1:length(res_eig)] <- res_eig
   return(eig_rand)
}


.permute_bootstrap_first_singval <- function(X, Y, length,
                                      from_indices_X, to_indices_X,
                                      from_indices_Y, to_indices_Y,
                                      center1, center2,
                                      scale1, scale2){
   # X and Y should be the original, non-preprocessed, non-MFA-normed data

   # First, get bootstrapped values for the first singular values of the subtables
   nN <- NROW(X)
   boot_index <- sample(nN, replace = TRUE)
   X_boot <- ExPosition::expo.scale(X[boot_index,], center1, scale1)
   Y_boot <- ExPosition::expo.scale(Y[boot_index,], center2, scale2)

   mfa_norm_res_X_boot <- purrr::map2(from_indices_X, to_indices_X, ~mfa_norm(X_boot, .x, .y))
   mfa_norm_res_Y_boot <- purrr::map2(from_indices_Y, to_indices_Y, ~mfa_norm(Y_boot, .x, .y))

   # bootstrapped first singular values
   singval_X <- purrr::map_dbl(mfa_norm_res_X_boot, 2)
   singval_Y <- purrr::map_dbl(mfa_norm_res_Y_boot, 2)

   # Now get a permuted sample
   ZX <- ExPosition::expo.scale(X, center1, scale1)
   ZY <- ExPosition::expo.scale(Y, center2, scale2)

   X_rand <- apply(ZX, 2, sample)
   Y_rand <- apply(ZY, 2, sample)

   # Normalize the random data by the bootstrapped first singular values
   mfa_norm_X_rand <- purrr::pmap(list(from_indices_X, to_indices_X, singval_X),
                                  ~ X_rand[,..1:..2] / ..3) %>%
      purrr::reduce(cbind)
   mfa_norm_Y_rand <- purrr::pmap(list(from_indices_Y, to_indices_Y, singval_Y),
                                  ~ Y_rand[,..1:..2] / ..3) %>%
      purrr::reduce(cbind)

   R_rand <- t(mfa_norm_X_rand) %*% mfa_norm_Y_rand

   eig_rand <- rep(0, length)
   if(dim(X)[2] < dim(Y)[2]){
      res_eig <- eigen(R_rand %*% t(R_rand), symmetric = TRUE, only.values = TRUE)$values
   } else {
      res_eig <- eigen(t(R_rand) %*% R_rand, symmetric = TRUE, only.values = TRUE)$values
   }

   eig_rand[1:length(res_eig)] <- res_eig
   return(eig_rand)

}
