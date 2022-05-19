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
add_knots=function(f,f_v = NULL,knots,L,M=5)
{
  # check the class of the data
  if (is.vector(f)){
    # transpose the vector to become a matrix
    f = matrix(data = f, nrow = 1, ncol = length(f))
  }

  nx=dim(f)[2] #The number of points in the grid
  #Evaluating AMSE for the input knots.
  K=length(knots)


  AMSE=vector('numeric',K-1) #Here the values of the average mean squared errors will be kept
  APPRERR=vector('numeric',L+1) #Here will be kept the sequence of improved approximation errors
  #(in the terms of the average squared L2 norm: ||f1 - hat f1||_2^2+...+||fn - hat fn||_2^2)
  #by piecewise constant functions resulting from adding the knots

  LE=knots[1:(K-1)] #The open ended left ends of the intervals (so add one to have the closed end)
  RE=knots[2:K]     #The close ended right ends of the intervals

  splits=vector('numeric',K-1) #The new interval-wise optimal split-points
  AMSE1=splits #The left-hand side (with respect to corresponding 'splits') values of the average mean square error
  AMSE2=splits #The right-hand side (with respect to corresponding 'splits') values of the average mean square error



  #First run through all the intervals is to compute all interval-wise split and corresponding 'AMSE1' and 'AMSE2'
  for(k in 1:(K-1)) #the loop running through all the intervals at the current knots values
  {
    #print(k)
    #browser()
    ff=f[,(knots[k]+1):(knots[k+1]), drop=FALSE]
    AMSE[k]=amse(ff) #Here we keep the average mean squared errors for the input knots
    newsp=opt_split(ff,AMSE[k],M=M) #Finding optimal split with the given interval
    splits[k]=knots[k]+newsp[[1]]
    AMSE1[k]=newsp[[2]]
    AMSE2[k]=newsp[[3]]
    AMSE_v <- c(newsp[[2]],newsp[[3]])
  }

  cat("proposed splits is", splits, "\n")
  # AMSE_v <- split(f_v, splits)
  # cat("first AMSE_v is", AMSE_v, "\n")
  # len_v = c(length(na.omit(f_v[1,0:splits])),length(na.omit(f_v[1,(splits+1):nx])))
  # #TMSE_v = c()
  # #nx_v = dim(f_v)[2]
  # nx_v <- length(na.omit(f_v[1,]))
  # print(nx_v)
  # cat("first len_v is", len_v, "\n" )
  # TMSE_v = sum(len_v*AMSE_v)/nx_v
  # cat("the first TMSE_V", TMSE_v, "\n")

  if(!is.null(f_v)){
    nx_v <- length(na.omit(f_v[1,]))
    # cat("nx is", nx_v, "\n")
    len_v = nx_v
    TMSE_v = c()
    AMSE_v = c()
  }


  #The average approximation error of the functions for the input set of knots.
  APPRERR[1]=sum((RE-LE)*AMSE)/nx  #Adding 1 is because the between knots interval is LE[i],RE[i]
  #so that the number of points in this interval is

  #computed splits (potential new knots) and corresponding AMSE1's and AMSE2's
  #The full set knots, while updated in the loop below are kept in
  FLE=LE #The initial values of the left endpoints
  FRE=RE #The initial values of the right endpoints
  FAMS=AMSE

  Fspl=splits    #The splits and corresponding AMS1 and AMS2
  FAMS1=AMSE1
  FAMS2=AMSE2

  for(i in 1:L){
    #START of the loop

    if(all(is.na(Fspl))){  #Checking if finding knots can be continued due to the constraint on the number M
      #of points per in between knots intervals
      warning(paste0('There are only ', K+i-1 ,' knots. Reduce L or M.' ))

      break
    }
    else{
      opt2=add_split(f,FLE,FRE,FAMS,FAMS1,FAMS2,Fspl, M=M)
      l=opt2$locsp #location of the new knot in the intervals used for the computation, i.e. the knot is
      #in (L[l]+1):R[l] so that the new intervals are (L[l]+1):NR[l], (NL[l+1]+1):R[l]
      #where NR[l]=NL[l+1] is the new split (knot)
      NL=opt2$NLE
      NR=opt2$NRE
      AMS=opt2$NAMSE

      APPRERR[1+i]=sum((NR-NL)*AMS)/nx #The new average sum of the squared norms of errors.

      AMS1=opt2$NAMSE1
      AMS2=opt2$NAMSE2
      spl=opt2$nsplits


      #Updating the complete set of knots by the knew knot.
      FLE=append(FLE,NL[l+1],after=l)
      FRE=append(FRE,NR[l],after=l-1)

      #Updating the average MSE
      FAMS[l]=AMS[l]
      FAMS=append(FAMS,AMS[l+1],after=l)

      #Updating the optimal splits
      Fspl[l]=spl[l]
      Fspl=append(Fspl,spl[l+1],after=l)

      #And the corresponding left and right average MSE
      FAMS1[l]=AMS1[l]
      FAMS1=append(FAMS1,AMS1[l+1],after=l)
      FAMS2[l]=AMS2[l]
      FAMS2=append(FAMS2,AMS2[l+1],after=l)

      cat("proposed splits is", spl, "\n")

      print("printing the new knot")
      print(NL[l+1])
      # calculating the amse for the validation set
      # cat("nsplits is", spl,"\n")
      # cat("l is",l,"\n")

      if(!is.null(f_v)){
        # calculate the current f_v
        c_f_v = f_v[,NL[l]:NR[l+1]]
        # print(dim(c_f_v))
        # print(NL)
        # print(NR)
        s_v = split(c_f_v,NR[l]-NL[l])
        if(NA %in% s_v){
          #browser()
          next
      }
        # length of the current f_v
        l_v = c(length(na.omit(f_v[1,(NL[l]+1):NR[l]])),length(na.omit(f_v[1,(NL[l+1]+1):NR[l+1]])))
        AMSE_v[l] = s_v[1]
        AMSE_v = append(AMSE_v, s_v[2],after = l)
        len_v[l] = l_v[1]
        len_v = append(len_v, l_v[2],after = l)
        # cat("AMSE_v is", AMSE_v,"\n")
        # cat("len_v is", len_v, sum(len_v), "\n")
        #TMSE_v = sum(len_v*AMSE_v)/nx_v
        TMSE_v = append(TMSE_v, sum(len_v*AMSE_v)/nx_v)

        cat("current TMSE is", TMSE_v, "\n")
        # browser()
      }

    }
    #END of the loop.
  }
  Fknots=c(FLE,FRE[length(FRE)])
  if(is.null(f_v)){
    add_knots=list(Fknots=Fknots,FAMSE=FAMS,APPRERR=APPRERR)
  } else{
    add_knots=list(Fknots=Fknots,FAMSE=FAMS,APPRERR=APPRERR,TMSE_v=TMSE_v)
  }

  return(add_knots)
}
########

