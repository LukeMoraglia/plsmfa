#' Create plot of loadings for two dimensions
#'
#' Creates a factor map from \code{PTCA4CATA::createFactorMap} using \code{axis1} on the horizontal
#'  and \code{axis2} on the vertical.
#'
#' @param loadings The matrix of loadings to plot
#' @param axis1 (default 1) column of \code{loadings} to plot on horizontal axis
#' @param axis2 (default 2) column of \code{loadings} to plot on vertical axis
#' @param col_points (default "blueviolet") either a single color, or a vector of colors of length \code{NROW{loadings}}
#' @param constraints (default NULL) a list of constraints for the axes. See \code{\link[PTCA4CATA]{createFactorMap}}
#' @param arrows (default TRUE) Whether to draw arrows from the origin to the points
#' @param ... args passed to \code{\link[PTCA4CATA]{createFactorMap}}
#'
#' @return A list containing:
#'   \item{\code{plot}}{The final plot created by the function}
#'   \item{\code{base_map}}{The output of \code{createFactorMap}}
#'   \item{\code{arrows}}{The output of \code{\link[PTCA4CATA]{addArrows}}}
#' @author Luke Moraglia
#' @export
loading_plot <- function(loadings, axis1 = 1, axis2 = 2,
                          col_points = "blueviolet", constraints = NULL,
                          arrows = TRUE, ...) {


   base_map <- PTCA4CATA::createFactorMap(loadings, axis1, axis2,
                              col.points = col_points,
                              col.labels = col_points,
                              constraints = constraints,
                              ...)

   if(arrows){
      the_arrows <- PTCA4CATA::addArrows(loadings, axis1, axis2,
                           color = col_points)
      final_plot <- base_map$zeMap + the_arrows
   } else{
      final_plot <- base_map$zeMap
   }
   return_list <- list(plot = final_plot,
                       base_map = base_map)
   if(arrows){
      return_list$arrows <- the_arrows
   }

   invisible(return_list)
}


#' Create a bar plot of bootstrap ratios (or any numerical vector)
#'
#' @param boot_ratios matrix of data to plot
#' @param dimension (default 1) which column to plot
#' @param threshold (default 2) bars below threshold are gray
#' @param above_threshold (default FALSE) whether to plot only the values above \code{threshold}
#' @param color_bar (default NULL) a vector of colors for the bars, the same length as the data
#' @param ... args passed to \code{\link[PTCA4CATA]{PrettyBarPlot2}}
#'
#' @return The bar plot, output from \code{\link[PTCA4CATA]{PrettyBarPlot2}}
#' @author Luke Moraglia
#' @export
boot_ratio_barplot <- function(boot_ratios, dimension = 1,
                               threshold = 2,
                               above_threshold = FALSE,
                               color_bar = NULL,
                               ...){

   barplot <- PTCA4CATA::PrettyBarPlot2(boot_ratios[,dimension],
                                        threshold = threshold,
                                        color4bar = color_bar,
                                        signifOnly = above_threshold,
                                        ...)

   invisible(barplot)
}


#' Create a factor map of two latent variables from the same matrix
#'
#' @param res_plsmfa output of \code{\link{plsmfa}}
#' @param axis1 (default 1) the latent variable to plot on the horizontal axis
#' @param axis2 (default 2) the latent variable to plot on the vertical axis
#' @param XorY (default "X") Either "X" or "Y", showing which matrix the LVs should come from
#' @param design (default NULL) A design vector for the rows. If included, group means are plotted
#' @param col_obs (default "blueviolet") Either one color or a vector of colors with a length of N observations
#' @param col_group (default "blueviolet") Either one color or a vector of colors with a length of N groups
#' @param n_iter (default 100) Number of bootstrap iterations for computing confidence intervals around the group means
#' @param xlab (default \code{paste("Latent Variable", axis1)}) Horizontal axis label
#' @param ylab (default \code{paste("Latent Variable", axis2)}) Vertical axis label
#' @param ... args passed to \code{\link[PTCA4CATA]{createFactorMap}} for the lv_fact_map
#'
#' @return A list containing:
#'   \item{\code{plot}}{The final plot from the function}
#'   \item{\code{lv_fact_map}}{Output from \code{\link[PTCA4CATA]{createFactorMap}} for the observations}
#'   \item{\code{labels}}{The axis labels}
#'   \item{\code{group_means}}{The group means on the latent variables}
#'   \item{\code{lv_group_means_fact_map}}{Output from \code{\link[PTCA4CATA]{createFactorMap}} for the group means}
#'   \item{\code{CI_ellipses}}{The confidence intervals of the group means}
#' @author Luke Moraglia
#' @export
latent_variable_map <- function(res_plsmfa,
                    axis1 = 1,
                    axis2 = 2,
                    XorY = "X",
                    design = NULL,
                    col_obs = "blueviolet",
                    col_group = "blueviolet",
                    n_iter = 100,
                    xlab = paste("Latent Variable", axis1),
                    ylab = paste("Latent Variable", axis2),
                    ...){
   inference <- TRUE
   if(is.null(design)){
      inference <- FALSE
   }
   if(XorY == "X"){
      lv <- res_plsmfa$pls$lx[,c(axis1, axis2)]
   } else {
      lv <- res_plsmfa$pls$ly[,c(axis1, axis2)]
   }
   column_names <- paste0("Dimension ", c(axis1, axis2))
   colnames(lv) <- column_names
   a_points <- 0.9
   if (inference) {
      a_points <- 0.2
   }
   lv_plot <- PTCA4CATA::createFactorMap(lv,
                                         col.points = col_obs,
                                         col.labels = col_obs,
                                         alpha.points = a_points,
                                         display.labels = FALSE,
                                         ...)
   labels = ggplot2::labs(x = xlab, y = ylab)
   if(inference){
      lv_group <- PTCA4CATA::getMeans(lv, design)
      lv_group_boot <- PTCA4CATA::Boot4Mean(lv, design, n_iter)
      colnames(lv_group_boot$BootCube) <- column_names
      mean_plot <- PTCA4CATA::createFactorMap(lv_group,
                                              col.points = col_group[rownames(lv_group)],
                                              col.labels = col_group[rownames(lv_group)],
                                              cex = 4,
                                              pch = 17,
                                              alpha.points = 0.8)
      CI_plot <- PTCA4CATA::MakeCIEllipses(lv_group_boot$BootCube[ ,c(1:2), ],
                                           col = col_group[rownames(lv_group)],
                                           names.of.factors = column_names)
      plot <- lv_plot$zeMap_background + lv_plot$zeMap_dots +
         mean_plot$zeMap_dots + mean_plot$zeMap_text + CI_plot + labels
   } else{
      plot <- lv_plot$zeMap + labels
   }

   result <- list(plot = plot,
                  lv_fact_map = lv_plot,
                  labels = labels)

   if(inference){
      result$group_means <- lv_group
      result$lv_group_means_fact_map <- mean_plot
      result$CI_ellipses <- CI_plot
   }

   invisible(result)


}

#' Latent variable map for a single dimension, X vs Y
#'
#'
#' @param res_plsmfa output of \code{\link{plsmfa}}
#' @param lv_num (default 1) which dimension's latent variables to plot
#' @param design (default NULL) A design vector for the rows. If included, group means are plotted
#' @param col_obs (default "blueviolet") Either one color or a vector of colors with a length of N observations
#' @param col_group (default "blueviolet") Either one color or a vector of colors with a length of N groups
#' @param n_iter (default 100) Number of bootstrap iterations for computing confidence intervals around the group means
#' @param xlab (default \code{paste0("Latent Variable X", lv_num)}) Horizontal axis label
#' @param ylab (default \code{paste0("Latent Variable Y", lv_num)}) Vertical axis label
#' @param ... args passed to \code{\link[PTCA4CATA]{createFactorMap}} for the lv_fact_map
#'
#' @return A list containing:
#'   \item{\code{plot}}{The final plot from the function}
#'   \item{\code{lv_fact_map}}{Output from \code{\link[PTCA4CATA]{createFactorMap}} for the observations}
#'   \item{\code{group_means}}{The group means on the latent variables}
#'   \item{\code{lv_group_means_fact_map}}{Output from \code{\link[PTCA4CATA]{createFactorMap}} for the group means}
#'   \item{\code{CI_ellipses}}{The confidence intervals of the group means}
#' @author Luke Moraglia
#' @export
latent_variable_XY_map <- function(res_plsmfa,
                                    lv_num = 1,
                                    design = NULL,
                                    col_obs = "blueviolet",
                                    col_group = "blueviolet",
                                    n_iter = 100,
                                    xlab = paste0("Latent Variable X", lv_num),
                                    ylab = paste0("Latent Variable Y", lv_num),
                                    ...){
   inference <- TRUE
   if(is.null(design)){
      inference <- FALSE
   }
   lv <- cbind(res_plsmfa$pls$lx[, lv_num], res_plsmfa$pls$ly[, lv_num])
   colnames(lv) <- c(xlab, ylab)
   a_points <- 0.9
   if (inference) {
      a_points <- 0.2
   }
   lv_plot <- PTCA4CATA::createFactorMap(lv,
                                      col.points = col_obs,
                                      col.labels = col_obs,
                                      alpha.points = a_points,
                                      display.labels = FALSE,
                                      ...)
   if(inference){
      lv_group <- PTCA4CATA::getMeans(lv, design)
      lv_group_boot <- PTCA4CATA::Boot4Mean(lv, design, n_iter)
      colnames(lv_group_boot$BootCube) <- c(xlab, ylab)
      mean_plot <- PTCA4CATA::createFactorMap(lv_group,
                                              col.points = col_group[rownames(lv_group)],
                                              col.labels = col_group[rownames(lv_group)],
                                              cex = 4,
                                              pch = 17,
                                              alpha.points = 0.8)
      CI_plot <- PTCA4CATA::MakeCIEllipses(lv_group_boot$BootCube[ ,c(1:2), ],
                                           col = col_group[rownames(lv_group)],
                                           names.of.factors = c(xlab, ylab))
      plot <- lv_plot$zeMap_background + lv_plot$zeMap_dots +
         mean_plot$zeMap_dots + mean_plot$zeMap_text + CI_plot
   } else{
      plot <- lv_plot$zeMap
   }

   result <- list(plot = plot,
                  lv_fact_map = lv_plot)

   if(inference){
      result$group_means <- lv_group
      result$lv_group_means_fact_map <- mean_plot
      result$CI_ellipses <- CI_plot
   }

   invisible(result)


}



#' Plot partial latent variable group means
#'
#' @param res_plsmfa output of \code{\link{plsmfa}}
#' @param design A design vector for the rows
#' @param lv_plot the output of \code{\link{latent_variable_map}}
#' @param axis1 (default 1) the latent variable to plot on the horizontal axis
#' @param axis2 (default 2) the latent variable to plot on the vertical axis
#' @param XorY (default "X") Either "X" or "Y", showing which matrix the LVs should come from
#' @param col_items (default NULL)  Either one color or a vector of colors with a length of N groups
#' @param col_blocks (default NULL)  Either one color or a vector of colors with a length of N subtables
#' @param use_defaults (default TRUE) set to FALSE to use other \code{...} args
#' @param ... args passed to \code{\link[PTCA4CATA]{createPartialFactorScoresMap}}
#'
#' @return A list containing:
#'   \item{\code{plot}}{The final plot from the function}
#'   \item{\code{part_lv_means_plot}}{Output of \code{\link[PTCA4CATA]{createPartialFactorScoresMap}}}
#'   \item{\code{lv_plot}}{The \code{lv_plot} from input}
#' @author Luke Moraglia
#' @export
partial_lv_group_means_map <- function(res_plsmfa,
                                       design,
                                       lv_plot,
                                       axis1 = 1,
                                       axis2 = 2,
                                       XorY = "X",
                                       col_items = NULL,
                                       col_blocks = NULL,
                                       use_defaults = TRUE,
                                       ...){
   if(XorY == "X"){
      lv <- res_plsmfa$pls$lx[, c(axis1, axis2)]
      part_lv <- res_plsmfa$partial_lv_X[ , c(axis1, axis2), ]
   } else {
      lv <- res_plsmfa$pls$ly[,c(axis1, axis2)]
      part_lv <- res_plsmfa$partial_lv_Y[ , c(axis1, axis2), ]
   }
   colnames(lv) <- colnames(part_lv) <- paste0("Dimension ", c(axis1, axis2))

   part_lv_means <- get_partial_lv_group_means(part_lv, design)

   if(use_defaults){
      part_lv_means_plot <- PTCA4CATA::createPartialFactorScoresMap(lv_plot$group_means,
                                                                   partialFactorScores = part_lv_means,
                                                                   colors4Items = col_items,
                                                                   colors4Blocks = col_blocks,
                                                                   alpha.lines = 0.4, alpha.points = 0.4,
                                                                   size.points = 3,
                                                                   size.labels = 4, shape.points = 18,
                                                                   alpha.labels = 0.6)
   } else {
      part_lv_means_plot <- PTCA4CATA::createPartialFactorScoresMap(lv_plot$group_means,
                                                                    partialFactorScores = part_lv_means,
                                                                    colors4Items = col_items,
                                                                    colors4Blocks = col_blocks,
                                                                    ...)
   }

   plot <- lv_plot$lv_fact_map$zeMap +
      lv_plot$lv_group_means_fact_map$zeMap_dots +
      lv_plot$lv_group_means_fact_map$zeMap_text +
      part_lv_means_plot$mapColByItems + lv_plot$labels

   return_list <- list(plot = plot,
                       part_lv_means_plot = part_lv_means_plot,
                       lv_plot = lv_plot
                       )
   invisible(return_list)
}

#' Create a correlation circle plot for either X or Y
#'
#' @param res_plsmfa output of \code{\link{plsmfa}}
#' @param axis1 (default 1) the latent variable to plot on the horizontal axis
#' @param axis2 (default 2) the latent variable to plot on the vertical axis
#' @param XorY (default "X") Either "X" or "Y", showing which matrix the LVs should come from
#' @param col_vars (default "blueviolet") Either one color or a vector of colors
#'  with a length of N variables
#'
#' @return A list containing:
#'  \item{\code{plot}}{The final plot from the function}
#'  \item{\code{base_plot}}{Output from \code{\link[PTCA4CATA]{createFactorMap}}}
#'  \item{\code{arrows}}{The output of \code{\link[PTCA4CATA]{addArrows}}}
#'  \item{\code{circle}}{The output of \code{\link[PTCA4CATA]{addCircleOfCor}}}
#' @author Luke Moraglia
#' @export
correlation_circle <- function(res_plsmfa,
                               axis1 = 1,
                               axis2 = 2,
                               XorY = "X",
                               col_vars = "blueviolet"){
   if(XorY == "X"){
      lv <- res_plsmfa$pls$lx[, c(axis1, axis2)]
      mat <- res_plsmfa$normed_X
   } else {
      lv <- res_plsmfa$pls$ly[,c(axis1, axis2)]
      mat <- res_plsmfa$normed_Y
   }

   colnames(lv) <- paste0("Dimension", c(axis1, axis2))
   cor <- stats::cor(mat, lv)
   cor_map <- PTCA4CATA::createFactorMap(cor,
                                         constraints = list(minx = -1,
                                                            maxx = 1,
                                                            miny = -1,
                                                            maxy = 1),
                                         col.labels = col_vars,
                                         display.points = FALSE,
                                         )
   the_arrows <- PTCA4CATA::addArrows(cor,
                                      color = col_vars)
   the_circle <- PTCA4CATA::addCircleOfCor()

   plot <- cor_map$zeMap + the_arrows + the_circle

   invisible(list(plot = plot,
               base_plot = cor_map,
               arrows = the_arrows,
               circle = the_circle))

}
