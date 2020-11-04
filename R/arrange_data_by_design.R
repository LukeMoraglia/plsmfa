#' Arrange columns of \code{data} into groups
#'
#' Rearranges the columns of a matrix based on \code{column_design}
#' so that the groups of columns are together
#'
#'
#' @param data data to rearrange
#' @param column_design a vector indicating group membership of columns
#'
#' @return a list with the data and column design rearranged into groups
#' @author Luke Moraglia
#' @export
arrange_data_by_design <- function(data, column_design){
   column_design <- as.matrix(column_design)
   design_mat <- ExPosition::makeNominalData(column_design)
   new_data <- list()
   new_column_design <- NULL
   for(i in 1:NCOL(design_mat)){
      new_data[[i]] = data[, design_mat[,i] == 1]
      new_column_design <- c(new_column_design,
                             column_design[design_mat[,i] == 1])
   }
   new_data <- do.call(cbind, new_data)
   return(list(data = new_data,
               column_design = new_column_design))
}
