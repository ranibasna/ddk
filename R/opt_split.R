#' @title Finding the optimal splitting point which minimize the average mean square error over the range of the data.
#' @description For a given data \code{f} and its average mean square error \code{AMSE}, the function returns a vector of
#'the optimal split with the average mean square errors on the left and right of the splitting point, respectively.
#' @param f  n x nx matrix. where n is number of samples (rows) and nx reperesents the grid size
#' in functional data interpreptation, AMSE is the total of mean square errors computed rowwise
#' @param AMSE  integer.  the average mean square errors of \code{f} computed over its range.
#' @param M integer.  the minimal number of points between the optimal split  and the two ends of the iterval.
#' The default is 10.  The program will return \code{Na} if there are less than \code{2*M+1} points in the range \code{f}.
#' @return A list made of three numeric values:  \code{opt_ix}, \code{AMSE_L} and \code{AMSE_R}.
#' The first numiric value \code{opt_ix} is the optimal split (knot), and  \code{AMSE_L}, \code{AMSE_R}
#' are the average mean square errors left and right the spliting point, respectively.
#' @export
#' @section References: Nassar, H., Podgórski, K. (2019) Empirically driven
#' orthonormal bases for functional data analysis. \emph{Preprint}. Department of Statistics, Lund University.
#' @examples ## Example:
#' n=10
#' f=rbetafda(n)
#' nx=dim(f)[2]
#' AMSE = amse(f) ### Total mean square error for whole the samples
#' opt_split(f,AMSE)
#' @seealso \code{split} for constructing split at a given knot;
#'

##################################

opt_split=function(f,AMSE, M=5){
  if (!(is.data.frame(f) | is.matrix(f) | is.vector(f))){
    stop("The argument x is not a dataframe, matrix, or vector")
  }
  if (!(is.numeric(AMSE) | is.integer(AMSE))){
    stop("The argument x is not a numeric or integer")
  }
  nx=dim(f)[2]   ### number of grid
  ## adding condition to avoid al NA inside one intervals moved to mse function
  if(nx<2*M+1){  # condition for having at least M points between the  optimal split  and the two ends of the iterval.
    opt_ix=NA
    AM1=NA
    AM2=NA
  }
  else
  {
    AM <- lapply(1:(nx-1),split,f=f)  # calculate the left and right AMSE for each discritization
    AM1 = unlist(lapply(AM, `[[`, 1)) # getting  the left  AMSE for each discritization
    AM2 = unlist(lapply(AM, `[[`, 2)) # getting  the right  AMSE for each discritization
    I = 1:(nx-1)
    DAMSE=AMSE* nx -(AM1*I+AM2*((nx-I)))
    DAMSE_new=DAMSE[(M+1):(nx-M-1)]
    if (all(DAMSE_new==0)){
      opt_ix= NA
      AM1=NA
      AM2=NA
    }
    else
    {
      ix=order(DAMSE_new, decreasing=TRUE) + M
      opt_ix= ix[1]
    }
  }
  opt_split=list(opt_ix,AM1[opt_ix],AM2[opt_ix])

  # names(opt_split)=c("optimal split","AMSE_L","AMSE_R")
  return(opt_split)
}
