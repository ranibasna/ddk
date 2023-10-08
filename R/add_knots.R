#' @title Add Quasi-optimal Knots for the Piecewise Constant Basis
#'
#' @description  Given discretized functional data \code{f} defined over a grid, find the knots that are quasi-optimal for
#' the piecewise constant basis. It is 'quasi-optimal' since the optimality is locally verified at each
#' step when a knot is added to already existing ones. There is no guarantee that the final set of knots
#' is optimal globally over all possible knots of this size.
#'
#'
#' Utilizing 'add_split' function for evaluation of L new knots only in the intervals that contain
#' more than the specified (by \code{M}) number of knots. If it is not possible for the function to return the specified by L
#' number of additional knots a warning message is reported informing what is the maximal number of new knots
#' that can be added. The knots are added in such a way that it achieves the largest drop in the AMSE (which the sum
#' of mean squared errors accross all functional data.)
#

#' @param f functional data, i.e. matrix of the values row-wise evaluated over equidistant arguments (not given)
#' @param f_v functional data i.e. matrix of the values row-wise evaluated over equidistant arguments as the validation data for which the validate the
#' decreas in the amse over the selected knots outputed from the functional data \code{f}
#' @param knots input knots, a sequence of ordered K integers in the range 0:nx, where nx=dim(f)[2]
#'              in between knots intervals are defined open ended on the left hand side (knots[i]) and closed
#'              on the right hand side (knots[i+1]): knots[i]+1,...,knots[i+1]
#' @param L number of additional knots
#' @param M the minimal number of points per intervals between knots for further split of the optimization to be performed.
#'          The default is 5. It means that if there is less than M points per interval the further split of the interval
#'          marked by knots is not performed. The program will stop if there is too few points to find the requested number
#'          of knots with the given restriction for M.
#'
#' @return A list of the following values
#'  \item{nknots}{sequence of K+L (old and new) knots}
#'  \item{NAMSE}{the corresponding sequence of K+L-1 of the within knots
#'                  average mean square errors}
#'  \item{APPRERR}{the decreasing sequence of the averaged squared L2 norms:
#'                 ||f1 - hat f1l||_2^2+...+||fn - hat fnl||_2^2, l=0,...,L, where hat f_il are piecewise constant
#'                 approximation of fi's with l knots added to the input knots.}
#' @export
#' @section References:
#' Nassar, H., Podgórski, K. Empirically driven orthonormal bases for functional data analysis.
#' \emph{Proceedings of European Numerical Mathematics and Advanced Applications Conference 2019}.
#' Cognitive Systems, Department of Applied Mathematics and Computer Science, Technical University of Denmark, Denmark
#'
#' Basna, R. Nassar, H., Podgórski, K. Machine Learning Assisted Orthonormal Bases Selection for Functional Data Analysis. (preprint)
#'
#' @examples
#' n=10 #number of samples
#' #f=rbetafda(n) #generating data
#' #f=rbetafda(n,ta=3,tb=3) #generating data
#' nx=1000
#' f=rbetafda(n,nx,ta=3,tb=3) #generating data
#'
#' nx=dim(f)[2] #size of the equidistant one dimensional grid
#' hf=1/(nx+1)  #increment s i z e
#' grid=matrix( seq (hf , 1-hf , by=hf) , nrow=1) #grid
#' xx=vector()
#' for( i in 1:(nx-1)){
#'   Q=split(f,i)
#'   xx[i]=Q[1]*i/nx+Q[2]*(nx-i)/nx
#' }
#' plot(xx)
#'
#' AMSE=c(mean((nx-1)/nx*apply(f,1,var)))
#' knots=c(0,nx) #We take zero as the location of the first knot since, we want intervals pointed
#' # by 'knots' to be open-close, i.e. the k-th interval is 'knots[k]+1, knots[k+1]'
#' K=length(knots)
#'
#' KS=add_knots(f,knots=knots,L=10, M = 5)
#' KS
#'
#' plot(log(KS$APPRERR))
#' @seealso \code{split}  for constructing split at a given knot; \code{opt_split}
#' for finding the optimal split within one interval; \code{add_splitw}
#' for selecting  the optimal split from a set of potential splits.
#'
#' @export


#############
add_knots = function(f, f_v=NULL, knots, L, M=5, auto_stop=FALSE, threshold=NULL, stop_method="absolute") {
  # check the class of the data
  if (is.vector(f)){
    # transpose the vector to become a matrix
    f = matrix(data = f, nrow = 1, ncol = length(f))
  }

  # Get the number of grid points
  nx = dim(f)[2]

  # Get the number of initial knots
  K = length(knots)


  # Initialize vectors for AMSE and approximation errors
  AMSE = vector('numeric', K-1)
  APPRERR = vector('numeric', L+1)
  # Get the left and right endpoints of the intervals
  LE = knots[1:(K-1)]
  RE = knots[2:K]
  # Initialize vectors for optimal split points and their corresponding AMSEs
  splits = vector('numeric', K-1)
  AMSE1 = splits
  AMSE2 = splits

  # Calculate initial AMSE and optimal split points for each interval
  for(k in 1:(K-1)) {
    ff = f[,(knots[k]+1):(knots[k+1]), drop=FALSE]
    AMSE[k] = amse(ff)
    newsp = opt_split(ff, AMSE[k], M=M)
    splits[k] = knots[k] + newsp[[1]]
    AMSE1[k] = newsp[[2]]
    AMSE2[k] = newsp[[3]]
  }
  # If validation set is provided, initialize validation-related variables
  cat("proposed splits is", splits, "\n")
  if(!is.null(f_v)){
    nx_v = length(na.omit(f_v[1,]))
    len_v = nx_v
    TMSE_v = c()
    AMSE_v = c()
  }

  # Calculate initial approximation error
  APPRERR[1] = sum((RE-LE) * AMSE) / nx

  # Copy initial knots and AMSE values
  FLE = LE
  FRE = RE
  FAMS = AMSE

  Fspl = splits
  FAMS1 = AMSE1
  FAMS2 = AMSE2

  prev_AMSE_v = NULL
  # Loop to add knots
  for(i in 1:L) {
    if(all(is.na(Fspl))) { #Checking if finding knots can be continued due to the constraint on the number M
      #of points per in between knots intervals
      warning(paste0('There are only ', K + i - 1, ' knots. Reduce L or M.'))
      break
    } else {
      opt2 = add_split(f, FLE, FRE, FAMS, FAMS1, FAMS2, Fspl, M=M)
      l = opt2$locsp #location of the new knot in the intervals used for the computation.
      NL = opt2$NLE
      NR = opt2$NRE
      AMS = opt2$NAMSE

      APPRERR[1+i] = sum((NR-NL) * AMS) / nx #The new average sum of the squared norms of errors.

      AMS1 = opt2$NAMSE1
      AMS2 = opt2$NAMSE2
      spl = opt2$nsplits

      #Updating the complete set of knots by the knew knot.
      FLE = append(FLE, NL[l+1], after=l)
      FRE = append(FRE, NR[l], after=l-1)

      #Updating the average MSE
      FAMS[l] = AMS[l]
      FAMS = append(FAMS, AMS[l+1], after=l)

      #Updating the optimal splits
      Fspl[l] = spl[l]
      Fspl = append(Fspl, spl[l+1], after=l)

      #And the corresponding left and right average MSE
      FAMS1[l] = AMS1[l]
      FAMS1 = append(FAMS1, AMS1[l+1], after=l)
      FAMS2[l] = AMS2[l]
      FAMS2 = append(FAMS2, AMS2[l+1], after=l)

      cat("proposed splits is", spl, "\n")

      print("printing the new knot")
      print(NL[l+1])
      # calculating the amse for the validation set
      if(!is.null(f_v) && auto_stop) {
        # calculate the current f_v
        c_f_v = f_v[, NL[l]:NR[l+1]]
        s_v = split(c_f_v, NR[l]-NL[l])

        if(NA %in% s_v) {
          next
        }
        # length of the current f_v
        l_v = c(length(na.omit(f_v[1, (NL[l]+1):NR[l]])), length(na.omit(f_v[1, (NL[l+1]+1):NR[l+1]])))
        AMSE_v[l] = s_v[1]
        AMSE_v = append(AMSE_v, s_v[2], after=l)
        len_v[l] = l_v[1]
        len_v = append(len_v, l_v[2], after=l)
        TMSE_v = append(TMSE_v, sum(len_v * AMSE_v) / nx_v)

        cat("current TMSE is", TMSE_v, "\n")
        if(!is.null(prev_AMSE_v)) {
          abs_diff = abs(TMSE_v[i] - prev_AMSE_v)

          if(stop_method == "absolute" && abs_diff < threshold) {
            break
          } else if(stop_method == "relative" && abs_diff < threshold * abs(TMSE_v[i])) {
            break
          }
        }

        prev_AMSE_v = TMSE_v[i]
      }
    }
  }

  Fknots = c(FLE, FRE[length(FRE)])
  if(is.null(f_v)){
    add_knots = list(Fknots=Fknots, FAMSE=FAMS, APPRERR=APPRERR)
  } else {
    add_knots = list(Fknots=Fknots, FAMSE=FAMS, APPRERR=APPRERR, TMSE_v=TMSE_v)
  }

  return(add_knots)
}
########
