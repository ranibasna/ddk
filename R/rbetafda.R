#' @title  Function to produce generic functional data
#' @param n number of samples
#' @param nx number of coloumns
#' @param ta is shape parameter
#' @param tb is shape parameter
#' @param seed is a seeding option for reproducibility issues
#' @return The function produces a fda random sample of size n over the equidistant grid of size nx, 1/(nx+1):1/(nx+1):nx/(nx+1)
#' with equal to square roots of beta densities with parameters alpha and beta simulated independently from gamma distributions with
#' the shape parameters ta and tb, respectively, and normalized to yield the mean 0.5.
#' @export
#' @importFrom  stats dbeta rgamma rnorm
#' @import graphics
#' @importFrom graphics grid
#'
#' @examples
#' f1=rbetafda(10)
#'

rbetafda=function(n,nx=500,ta=1,tb=1,sigma=0,K=1000, seed)
{
  if(missing(seed)){
    set.seed(11)
  }else{
    set.seed(seed)
  }

  grid=matrix(seq(1/(nx+1),nx/(nx+1),by=1/(nx+1)),nrow=1)
  f=matrix(0,ncol=nx,nrow=n)
  BB=matrix(0,ncol=nx,nrow=1)
  for(j in 1:n)
  {
    if(sigma>0){
      ZZ=as.matrix(rnorm(K),ncol=1)
      KK=matrix(1:K,ncol=1)
      lambda=1/KK^2/pi^2
      E=sqrt(2)*sin(pi*KK%*%grid) #a row contains values of an eigenfunction
      BB=matrix(0,ncol=nx,nrow=1)
      for(i in 1:K)
      {
        BB[1,]=sqrt(lambda[i])*E[i,]*ZZ[i]+BB[1,]
      }
    }
    a=rgamma(1,ta,1/ta)/2
    b=rgamma(1,tb,1/tb)/2
    f[j,]=sqrt(dbeta(t(grid),a,b))+BB[1,]
  }
  f
}


# f1=rbetafda(10)
#
#
# nx=dim(f1)[2] #size of the equidistant one dimensional grid
# hf=1/(nx+1)  #increment s i z e
# grid=matrix( seq (hf , 1-hf , by=hf) , nrow=1) #grid
#
# par(mfrow=c(2,2))
# plot(grid,f1[1,],type='l',ylab='',xlab='',ylim=c(min(f1),max(f1)))
# for(i in 2:10)
# {
# lines(grid,f1[i,],type='l',col=i)
# }


