#' @title adding an optimal splitting point by choosing it from a set of potential splits.
#'
#' @description for a given data \code{f}, and a set of intervals (that have end points given in 'LE' and 'RE'),
#'  left and right hand side endpoints,
#' all in the range [1:nx], this function finds the an additional split (knot) from the set of potential knots \code{splits}
#'  by minimizing the average mean square error of piecewise constant approximation to (rowwise) functional data given in f
#' that result from adding a new split. The intervals do not need to cover the entire range of f.
#' We use convention that for LE[k], the k-th interval is '(LE[k]+1):RE[k]'
#' which corresponds to open-closed intervals in the 'continuous' interpretation of the intervals (LE[k],RE[k]]
#' The minimization is over all k. The noimnated optimal splits within each interval are already given in the input and
#' the algorithm returns new optimal within interval splits that are increased by one and only modified in
#' the two intervals
#'
#'
#' @param f  n x nx matrix, where n is number of samples (rows) and nx reperesents the grid size in functional data interpreptation.
#' @param M the minimal number of points per intervals between knots for further split of the optimization to be performed.
#'          The default is 5. It means that if there is less than M points per interval the further split of the interval
#'          marked by knots is not performed. The program will stop if there is too few points to find the requested number
#'          of knots with the given restriction for M.
#' @param LE vector of integers. indices for the left endpoints that are given as vectors of K integers in the range from 0 to nx-1,
#' @param RE vector of integers. indices for the right endpoints that are given as vectors of K integers in the range from 1 to nx,
# They are ordered and moreover: LE[k]<RE[k]<=LE[k+1]<RE[k+1], k in 1:K (do not need to cover the whole range of f)
#' @param AMSE vector of K values of the average mean square errors of the piecewise constant approximation based on given knots
#' @param AMSE1 vector of K values of the left-hand side (with respect to corresponding 'splits' below) values of the average mean square error
#' @param AMSE2 vector of K values of the righ-hand side (with respect to corresponding 'splits' below) values of the average mean square error
#' @param splits vector of K values of the suggested splits (they are chosen optimally within given interval and
#' will be used to choose a new knot). The knot is chosen according to the criterium: the split leading to the largest gain accross all between knots intervals becomes a new knot.
#'
#' @return the list of the following:
#' \itemize{
#' \item{NLE} and \item{NRE} -- two vectors of all K+1 ordered indices within the range [1:nx]representing end points of new increased by one set of intervals obtained
#' by adding the optimal knot (the optimal split of one of the input intervals)
#' \item{nsplits}-{a vector of K+1 values of the suggested splits (for using in the next iteration of this functiom to choose a new knot). These splits are optimal within intervals.}
#' \item{NAMSE}-{the K+1 vector of new the average mean square errors over the intervals obtained from the knots including the new one.}
#' \item{NAMSE1}-{the vector of K+1 values of the left-hand side (with respect to corresponding 'splits') values of the average mean square errors}
#' \item{NAMSE2}-{a vector of K+1 values of the right-hand side (with respect to corresponding 'splits') values of the average mean square errors}
#' \item{locsp}-{the location (left hand side) of the added new optimal split location,i.e. the split occurs in (LE[k]+1):RE[k], where k=locsp while LE and RE are the input values,
#' in the terms of output the new split is at NLE[k+1]=NRE[k]}
#' }
#' @export
#' @section References: Nassar, H., Podgórski, K. (2019) Empirically driven
#'   orthonormal bases for functional data analysis. \emph{Preprint}. Department
#'   of Statistics, Lund University.
#'
#' @examples
#'#START of example
#' #Generating data, needs 'rbetafda' function
#' n=10
#' f=rbetafda(n)
#' nx=dim(f)[2] #size of the equidistant one dimensional grid
#' hf=1/(nx+1)  #increment s i z e
#' grid=matrix( seq (hf , 1-hf , by=hf) , nrow=1) #grid
#' xx=vector()
#' for( i in 1:499){
#'  Q=split(f,i)
#'  xx[i]=Q[1]*i/nx+Q[2]*(nx-i)/nx
#'  }
#'  plot(xx)
#'  #Preparation of the input values
#' AMSE=c(mean((nx-1)/nx*apply(f,1,var)))
#' knots=c(0,nx) #We take zero as the location of the first knot since, we want
#' # intervals pointed by 'knots' to be open-close, i.e. the k-th interval is 'knots[k]+1, knots[k+1]'
#' K=length(knots)
#' LE=knots[1:(K-1)] #The open ended left ends of the intervals (so add one to have the closed end)
#' RE=knots[2:K]     #The close ended right ends of the intervals
#' splits=vector('numeric',K-1) #The new interval-wise optimal split-points
#' #The left-hand side (with respect to corresponding 'splits') values of the average mean square error
#' AMSE1=splits
#' # The right-hand side (with respect to corresponding 'splits')
#' # values of the average mean square error
#' AMSE2=splits
#' # First run through all the intervals is to compute all
#' # interval-wise split and corresponding 'AMSE1' and 'AMSE2'
#' for(k in 1:(K-1)) #the loop running through all the intervals at the current knots values
#' {
#'   ff=f[,(knots[k]+1):(knots[k+1])]
#'   newsp=opt_split(ff,AMSE[k])
#'   splits[k]=knots[k]+newsp[[1]]
#'   AMSE1[k]=newsp[[2]]
#'   AMSE2[k]=newsp[[3]]
#'  }
#' # computed splits (potential new knots) and corresponding AMSE1's and AMSE2's
#' ####################
#' opt2=add_split(f,LE,RE,AMSE,AMSE1, AMSE2, splits)
#'  L=opt2$NLE
#'  R=opt2$NRE
#'  AMS=opt2$NAMSE
#'  AMS1=opt2$NAMSE1
#'  AMS2=opt2$NAMSE2
#'  spl=opt2$nsplits
#'  opt2=add_split(f,L,R,AMS,AMS1,AMS2,spl)


add_split=function(f,LE,RE,AMSE,AMSE1,AMSE2,splits, M=5)
{
  K=length(AMSE) #The common size of all vectors it will increase by one at the end of the function

  #Next, one of the splits (optimal) becomes a new knot (and thus one endpoint is added to LE and RE each).

  nLE=LE  #new endpoints they will increase by one in the next part
  nRE=RE  #of the code
  NAMSE=AMSE #New AMSE corresponding to 'nknots'

  #Compute the drops in the average mean square errors, order them, and get ranks.

  #DAMSE=AMSE-((splits-LE)/(RE-LE))*AMSE1-((RE-splits)/(RE-LE))*AMSE2
  DAMSE=AMSE*(RE-LE)-(splits-LE)*AMSE1-(RE-splits)*AMSE2
  #DAMSE=((splits-LE)/(RE-LE))*AMSE1+((RE-splits)/(RE-LE))*AMSE2
  cat("DMSE add_split", DAMSE, "\n")


  pos=which.max(DAMSE) #Location of the optimal split (also left hand side endpoint of the corresponding knot-interval)
  nk=splits[pos] #The split leading to the largest gain, this becomes a new knot.

  nLE=append(nLE,nk,after=pos) #Now we have K+1 knots, so one more left hand side endpoint
  nRE=append(nRE,nk,after=(pos-1)) # and one more right hand side endpoint


  #After getting a new knot we need change also AMSE
  NAMSE[pos]=AMSE1[pos]
  NAMSE=append(NAMSE,AMSE2[pos],after=pos) #Now we have K AMSE's

  #Now it is time to update splits and AMSE1 and AMSE2
  nsplits=splits #New splits corresponding to 'nknots' (to be revised below)
  NAMSE1=AMSE1 #New AMSE1 corresponding to 'nknots' (to be revised below)
  NAMSE2=AMSE2 #New AMSE1 corresponding to 'nknots' (to be revised below)

  #Compute the optimal splits on the left and right intervals around the new knot
  ff1=f[,(nLE[pos]+1):nRE[pos], drop=FALSE]  #Part of the data below the new knot, the new knot is at LE[pos+1]=RE[pos]
  newsp=opt_split(ff1,NAMSE[pos], M = M)

  nsplits[pos]=nLE[pos]+newsp[[1]]       #Computing the split, AMSE1, and AMSE2 for the split below the new knot
  NAMSE1[pos]=newsp[[2]]
  NAMSE2[pos]=newsp[[3]]

  ff2=f[,(nLE[pos+1]+1):nRE[pos+1], drop=FALSE] #Computing the split, AMSE1, and AMSE2 for the split above the new knot
  newsp=opt_split(ff2,NAMSE[pos+1], M = M)

  nsplits=append(nsplits,nLE[pos+1]+newsp[[1]],after=pos) #Computing the split, AMSE1, and AMSE2 for the split above the new knot
  NAMSE1=append(NAMSE1,newsp[[2]],after=pos)
  NAMSE2=append(NAMSE2,newsp[[3]],after=pos)

  opt_knots=list(NLE=nLE,NRE=nRE,nsplits=nsplits,NAMSE=NAMSE,NAMSE1=NAMSE1,NAMSE2=NAMSE2,locsp=pos)
  return(opt_knots)
}


