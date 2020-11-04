#' Bootstrap for loadings of PLSMFA
#'
#' \code{boot_plsmfa} performs bootstrap resampling to generate a sampling distribution of the loadings.
#' This distribution can be used to generate confidence intervals, or bootstrap ratios.
#'
#' Bootstrap resampling is used to assess the stability of the results by simulating how the results would
#' vary if the data were collected again from the same population as the sample data. It simulates this process
#' by sampling the observations of the data \emph{with replacement}. This "bootstrapped sample" is analyzed
#' with PLSMFA and the results of the loadings are recorded. \code{boot_plsmfa} repeats these steps
#' \code{n_iter} times and stores the "cube" of loadings. From this cube, confidence intervals could be
#' calculated, but more
#' often (especially when there is a large number of variables), a single value called a bootstrap ratio is
#' used to assess the stability of the loading. The bootstrap ratio is similar to a $t$-value, since it is
#' the mean of the sampling distribution divided by the standard deviation of the distribution.
#'
#' At present, \code{boot_plsmfa} computes bootstrap ratios, but not confidence intervals. Also, as a technical
#' aside, because the SVD is prone to nuissance variance from reflection, reordering, and rotation (see Milan
#' and Whittaker, 1995), \code{boot_plsmfa} implements a Procrustes rotation on the loadings of each iteration
#' (\code{uv_boot}).
#'
#' @inheritParams plsmfa
#' @param n_iter (default 1000) The number of iterations for the bootstrap
#' @param n_dimensions (default 0) The number of dimensions to bootstrap. The default 0 means to use all
#'  of the dimensions. For data with many variables, using all of the dimensions will take longer to compute.
#' @param boot_ratio_threshold (default 2) The threshold for deeming bootstrap ratios as "significant"
#'
#' @return A list that contains:
#'   \item{\code{u_boot_cube}}{A 3D array containing \code{n_iter} bootstrapped versions of \code{u}}
#'   \item{\code{boot_ratios_i}}{Bootstrap ratios for the \eqn{I} variables of \code{data1}}
#'   \item{\code{sig_boot_ratios_i}}{Whether each bootstrap ratio for the \eqn{I} variables of \code{data1}
#'    is above the threshold.}
#'   \item{\code{v_boot_cube}}{A 3D array containing n_iter bootstrapped versions of \code{v}}
#'   \item{\code{boot_ratios_j}}{Bootstrap ratios for the J variables of \code{data2}}
#'   \item{\code{sig_boot_ratios_j}}{Whether each bootstrap ratio for the J variables of \code{data2}}
#'
#'
#' @author Luke Moraglia
#' @export
#'
#'
boot_plsmfa <- function(data1, data2,
                        column_design1, column_design2,
                        center1 = TRUE, center2 = TRUE,
                        scale1 = "SS1", scale2 = "SS1",
                        n_iter = 1000,
                        n_dimensions = 0,
                        boot_ratio_threshold = 2){
   # Should probably add more checks like these to all the main methods
   nN <- NROW(data1)
   if(nN != NROW(data2)){
      stop("data1 and data2 do not have the same number of rows")
   }

   nI <- NCOL(data1)
   nJ <- NCOL(data2)
   maxRank <- min(nI,nJ)
   if (maxRank < n_dimensions){
      n_dimensions <- maxRank
   }

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

   res_pls <- first_pos_gplssvd(GSVD::gplssvd(mfa_norm_X, mfa_norm_Y, k = n_dimensions), mfa_norm_X, mfa_norm_Y)

   if (n_dimensions > length(res_pls$d) | n_dimensions == 0){
      n_dimensions <- length(res_pls$d)
   }

   # I-set
   u_boot_cube <- array(NA, dim = c(nI, n_dimensions, n_iter))
   dimnames(u_boot_cube)[1] <- list(colnames(data1))
   dimnames(u_boot_cube)[2] <- list(paste0("Dimension ", 1:n_dimensions))
   dimnames(u_boot_cube)[3] <- list(paste0("Iteration ", 1:n_iter))

   # J-set
   v_boot_cube <- array(NA, dim = c(nJ, n_dimensions, n_iter))
   dimnames(v_boot_cube)[1] <- list(colnames(data2))
   dimnames(v_boot_cube)[2] <- list(paste0("Dimension ", 1:n_dimensions))
   dimnames(v_boot_cube)[3] <- list(paste0("Iteration ", 1:n_iter))

   for(i in 1:n_iter){
      boot.index <- sample(nN, replace = TRUE)
      X_boot <- ExPosition::expo.scale(data1[boot.index,], center1, scale1)
      Y_boot <- ExPosition::expo.scale(data2[boot.index,], center2, scale2)

      mfa_norm_res_X_boot <- purrr::map2(from_indices_X, to_indices_X, ~mfa_norm(X_boot, .x, .y))
      mfa_norm_res_Y_boot <- purrr::map2(from_indices_Y, to_indices_Y, ~mfa_norm(Y_boot, .x, .y))

      mfa_norm_X_boot <- purrr::map(mfa_norm_res_X_boot, 1) %>% purrr::reduce(cbind)
      mfa_norm_Y_boot <- purrr::map(mfa_norm_res_Y_boot, 1) %>% purrr::reduce(cbind)

      colnames(mfa_norm_X_boot) <- colnames(ZX)
      colnames(mfa_norm_Y_boot) <- colnames(ZY)

      # faster than calling gplssvd because we only need the loadings
      R_boot <- t(mfa_norm_X_boot) %*% mfa_norm_Y_boot
      res_svd_boot <- GSVD::tolerance_svd(R_boot, nu = n_dimensions, nv = n_dimensions)

      # puts the loadings into one matrix (I + J) by L
      uv_boot <- rbind(res_svd_boot$u, res_svd_boot$v)

      # Procrustes rotation, see Milan and Whittaker, 1995
      uv_orig <- rbind(res_pls$u, res_pls$v)
      nI <- NROW(res_pls$u)
      nJ <- NROW(res_pls$v)
      NOP <- GSVD::tolerance_svd(t(uv_boot) %*% uv_orig)
      q <- NOP$u %*% t(NOP$v)

      uv_boot_hat <- uv_boot %*% q

      u_boot_cube[,,i] <- uv_boot_hat[1:nI,]
      v_boot_cube[,,i] <- uv_boot_hat[(nI + 1):(nI + nJ),]
   }

   # compute bootstrap ratios
   BR_i <- .boot_ratio(u_boot_cube, boot_ratio_threshold)
   BR_j <- .boot_ratio(v_boot_cube, boot_ratio_threshold)

   result <- list(u_boot_cube = u_boot_cube,
                  boot_ratios_i = BR_i$boot_ratios,
                  sig_boot_ratios_i = BR_i$sig_boot_ratios,
                  v_boot_cube = v_boot_cube,
                  boot_ratios_j = BR_j$boot_ratios,
                  sig_boot_ratios_j = BR_j$sig_boot_ratios,
                  fixed_res = res_pls)

   invisible(result)

}



.boot_ratio <- function(boot_cube, boot_ratio_threshold = 2){
   # Code inspired by data4PCCAR::Boot4PLSC
   # Calculate mean of each loading sampling distribution
   boot_cube_mean <- apply(boot_cube, c(1,2), mean)
   # Calculate sd of each loading sampling distribution
   boot_cube_sd <- apply(boot_cube, c(1,2), stats::sd)
   boot_ratios <- boot_cube_mean / boot_cube_sd
   sig_boot_ratios <- (abs(boot_ratios) > boot_ratio_threshold)
   rownames(boot_ratios) <- rownames(boot_cube)
   rownames(sig_boot_ratios) <- rownames(boot_cube)
   return(list(boot_ratios = boot_ratios,
               sig_boot_ratios = sig_boot_ratios))
}

