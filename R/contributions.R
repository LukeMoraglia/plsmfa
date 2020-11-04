#' Compute contributions for variables and subtables
#'
#' The function uses the left singular vectors (loadings for X), the right singular
#' vectors (loadings for Y), and the column design vectors to compute the
#' contributions for each variable and each subtable.
#'
#' @param left_load left singular vectors / loadings for X
#' @param right_load right singular vectors / loadings for Y
#' @param column_design_X design vector for columns of X
#' @param column_design_Y design vector for columns of Y
#'
#' @return A list containing:
#'  \item{\code{ctr_var_X}}{Contributions for variables of X}
#'  \item{\code{ctr_var_Y}}{Contributions for variables of Y}
#'  \item{\code{ctr_sub_X}}{Contributions for subtables of X}
#'  \item{\code{ctr_sub_Y}}{Contributions for subtables of Y}
#'
#' @author Luke Moraglia
#' @export
contributions <- function(left_load, right_load, column_design_X, column_design_Y){
   # ctrs have variables on the rows and dimensions on the columns
   ctr_var_X <- left_load^2
   ctr_var_Y <- right_load^2
   # The getMeans function doesn't actually get the means. It gets the sum
   # for each column based on the groups of variables in the column design.
   # ctr_sub have rows corresponding to groups of variables and columns
   # with dimensions
   ctr_sub_X <- PTCA4CATA::getMeans(ctr_var_X, column_design_X, sum)
   ctr_sub_Y <- PTCA4CATA::getMeans(ctr_var_Y, column_design_Y, sum)
   return(list(ctr_var_X = ctr_var_X,
               ctr_var_Y = ctr_var_Y,
               ctr_sub_X = ctr_sub_X,
               ctr_sub_Y = ctr_sub_Y))

}
