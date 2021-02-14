#--- Functions are grouped with respect to the part of the package they are used ---#
#-------------------- Within the group the order is hierarchical  ---------------------#

#------------------------#
#'--- AuxFun for adding knots ---#
#' ------------------------#




#' Function mean squared error 'mse'
#' For a given vector x, the function 'mse' returns the the mean squared error,
#' mse measures the average of the squares of the errors
#'
#' @param x data of an n x nx matrix, where n is number of samples (rows) and nx reperesents the grid size
# in functional data interpreptation,
#' @export
#' @importFrom stats na.omit

mse=function(x){
  if (!(is.data.frame(x) | is.matrix(x) | is.vector(x))){
    stop("The argument x is not a dataframe, matrix or vector")
  }
  if(all(is.na(x))){
    return(NA)
  }
  else{
    x <- na.omit(x)
    return(sum((x-mean(x))^2)/length(x))
  }
}


#' Function average Mean Square Error 'AMSE'
#' For a given data f, the function 'AMSE' returns the
#' average mean square error for the data over all samples
#' @param f data of an n x nx matrix, where n is number of samples (rows) and nx reperesents the grid size
# in functional data interpreptation.
#' @export
#' @importFrom stats na.omit

amse= function(f){
  if (!(is.data.frame(f) | is.matrix(f))){
      stop("The argument x is not a dataframe, matrix")
 }
  # mse
  row_mse = apply(f,1,mse)
  if(all(is.na(row_mse))){
    return(NA)
  }
  else{
    return(mean(row_mse, na.rm = TRUE))
  }
}




