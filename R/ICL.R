#' Finds the active model's parameters
#'
#' @param M matrix
#' @param n number of units (equal to number of \code{M}'s rows)
#' @param k parameters to retrieve
#' @param na.rm logical, whether to ignore \code{NA} data
#'
#' @return An array containing the parameters
#'
#' @export
findActiveParam<-function(M, n, k, na.rm=TRUE){
  parByHB<-outer(k,k)
  parByHB<-array(parByHB,dim=c(dim(parByHB),dim(M)[3]))
  if(sum(n)!=dim(M)[1]) stop("Dimensions do not match!")
  kVec<-1:length(n)
  clu<-rep(kVec, n)
  for(r in 1:dim(M)[3]){
    for(i in kVec)for(j in kVec){
      parByHB[i, j,r]<-parByHB[i, j,r]*(var(as.vector(M[clu==i, clu==j,r]), na.rm = na.rm)>0)
    }
  }
  return(parByHB)
}



#' Function that computes integrated classification likelihood based on stochastic one-mode and linked block modeling. If \code{clu} is a list, the method for linked/multilevel networks is applied. The support for multirelational networks is not tested.
#'
#' @param M A matrix representing the (usually valued) network. For multi-relational networks, this should be an array with the third dimension representing the relation.
#' @param clu A partition. Each unique value represents one cluster. If the nework is one-mode, than this should be a vector, else a list of vectors, one for each mode. Similarly, if units are comprised of several sets, clu should be the list containing one vector for each set.
#' @param weights The weights for each cell in the matrix/array. A matrix or an array with the same dimmensions as \code{M}.
#' @param uWeights The weights for each unin. A vector with the length equal to the number of units (in all sets).
#' @param diagonal How should the diagonal values be treated. Possible values are:
#' \itemize{
#'   \item ignore - diagonal values are ignored
#'   \item seperate - diagonal values are treated seperately
#'   \item same - diagonal values are treated the same as all other values
#' }
#' @param limitType Type of limit to use. Forced to 'none' if \code{limits} is \code{NULL}. Otherwise, one of either \code{outer} or \code{inner}.
#' @param limits If \code{diagonal} is \code{"ignore"} or \code{"same"}, an array with dimensions equal to:
#' \itemize{
#'   \item number of clusters (of all types)
#'   \item number of clusters (of all types)
#'   \item number of relations
#'   \item 2 - the first is lower limit and the second is upper limit
#' }
#' If \code{diagonal} is \code{"seperate"}, a list of two array. The first should be as described above, representing limits for off diagonal values. The second should be similar with only 3 dimensions, as one of the first two must be omitted.
#' @param addOne Should one tie with the value of the tie equal to the density of the superBlock be added to each block to prevent block means equal to 0 or 1 and also "shrink" the block means toward the superBlock mean. Defaults to TRUE.
#' @param eps If addOne = FALSE, the minimal deviation from 0 or 1 that the block mean/density can take.
#' @param weightClusterSize The weight given to cluster sizes (logprobabilites) compared to ties in loglikelihood. Defaults to 1, which is "classical" stochastic blockmodeling.
#'
#' @return The value of ICL
#' 
#' @examples 
#' # Create a synthetic network matrix
#' set.seed(2022)
#' library(blockmodeling)
#' k<-2 # number of blocks to generate
#' blockSizes<-rep(20,k)
#' IM<-matrix(c(0.8,.4,0.2,0.8), nrow=2)
#' clu<-rep(1:k, times=blockSizes)
#' n<-length(clu)
#' M<-matrix(rbinom(n*n,1,IM[clu,clu]),ncol=n, nrow=n)
#' clu<-sample(1:2,nrow(M),replace=TRUE)
#' plotMat(M,clu) # Have a look at this random partition
#' ICL_pre<-ICLStochBlock(M,clu) # Calculate its ICL
#' ICL_pre
#' res<-stochBlock(M,clu=clu) # Optimizing the partition
#' plot(res) # Have a look at the optimized partition
#' ICL_post<-res$ICL # Calculate its ICL
#' ICL_post
#' # We expect the ICL pre-optimisation to be smaller:
#' ICL_pre<ICL_post			   
#' 
#' @seealso \code{\link{llStochBlock}}; \code{\link{weightsMlLoglik}}
#'
#' @export
ICLStochBlock<-function(M,
                       clu,
                       weights=NULL,
                       uWeights=NULL,
                       diagonal = c("ignore","seperate","same"),
                       limitType=c("none","inside","outside"),
                       limits=NULL,
                       weightClusterSize=1.0,
					   addOne = TRUE,
					   eps = 0.001){

  n1<-dim(M)[1]
  if(is.list(clu)) {
    n<-sapply(clu, length)
    k<-sapply(clu, function(x)length(unique(x)))
  }  else{
    n<-length(clu)
    k<-length(unique(clu))
  }

  if(sum(n)!=n1) stop("The length of clu and dimension of M does not match!")
  diagonal<-match.arg(diagonal)
  limitType<-match.arg(limitType)
  if(is.null(weights)){
    weights<-M
    weights[]<-1
  } else if(any(dim(weights)!=dim(M))) stop("Weights have wrong dim!")
  w<-weights
  if(is.null(uWeights)){
    uWeights<-rep(1.0, n1)
  }
  if(length(uWeights)!=n1) stop("uWeights has wrong length!")
  
  
  nMode<-ifelse(is.list(clu),length(clu),1)
  
  if(nMode>1){
    tmN<-sapply(clu,length)
    clu<-lapply(clu,function(x)as.integer(factor(x)))
    tmNclu<-sapply(clu,max)
    for(iMode in 2:nMode){
      clu[[iMode ]]<-clu[[iMode ]]+sum(tmNclu[1:(iMode -1)])
    }
    clu<-unlist(clu)
  } else{
    clu<-as.integer(factor(clu))
    tmNclu<-max(clu)
    tmN<-length(clu)
  }
  clu <- clu - 1
  if(length(dim(M))==2) M<-array(M,dim=c(dim(M),1))
  if(length(dim(w))==2) w<-array(w,dim=c(dim(w),1))
  
  if(is.null(limits)){
    bordersMatLower <- bordersMatUpper <- bordersSeperateLower <- bordersSeperateUpper<-NULL
    if(limitType!="none"){
      limitType<-"none"
      warning("limitType is set to 'none' as limits are NULL!")
    }
  } else {
    if(diagonal %in% c("ignore","same")){
      bordersSeperateLower <- bordersSeperateUpper
      if(is.list(limits)){
        limits<-limits[[1]]
      }
      if(!is.array(limits)) stop("'limits' must be specified as an array!")
      dl<-dim(limits)
      if(length(dl)==3 & dim(M)[3]==1) limits<-array(limits, dim=c(dl[1:2],1,dl[3]))
      if(dim(limits)!=4) stop("'limits' has wrong dimmensions (see help for correct dimmensions)")
      
      if(all(dim(limits)!=c(sum(tmNclu),sum(tmNclu),dim(M)[3],2))){
        stop("'limits' has wrong dimmensions (see help for correct dimmensions)")
      } else{
        bordersMatLower <- limits[,,,1]
        bordersMatUpper <- limits[,,,2]
      }
    } else {
      if(is.list(limits) & length(limits)==2){
        diagLimits<-limits[[2]]
        limits<-limits[[1]]
      } else stop("If diagonal is 'seperate', limits must be a list of length 2")
      if(!is.array(limits)) stop("First element of 'limits' must be specified as an array!")
      dl<-dim(limits)
      if(length(dl)==3 & dim(M)[3]==1) limits<-array(limits, dim=c(dl[1:2],1,dl[3]))
      if(dim(limits)!=4) stop("First element of 'limits' has wrong dimmensions (see help for correct dimmensions)")
      if(all(dim(limits)!=c(sum(tmNclu),sum(tmNclu),dim(M)[3],2))){
        stop("First element of 'limits' has wrong dimmensions (see help for correct dimmensions)")
      } else{
        bordersMatLower <- limits[,,,1]
        bordersMatUpper <- limits[,,,2]
      }
      
      if(!is.array(diagLimits)) stop("Second element of 'limits' must be specified as an array!")
      dl<-dim(diagLimits)
      if(length(dl)==2 & dim(M)[3]==1) limits<-array(limits, dim=c(dl[1],1,dl[2]))
      if(dim(diagLimits)!=3) stop("Second element of 'limits' has wrong dimmensions (see help for correct dimmensions)")
      if(all(dim(diagLimits)!=c(sum(tmNclu),dim(M)[3],2))){
        stop("Second element of 'limits' has wrong dimmensions (see help for correct dimmensions)")
      } else{
        bordersSeperateLower <- diagLimits[,,,1]
        bordersSeperateUpper <- diagLimits[,,,2]
      }
    }
  }
  
  if(diagonal == "ignore")for(i in 1:dim(w)[3]){
    diag(w[,,i])<-0
  }
  w<-w*findEmptySuperbocks(M,n = n)
  w<-w/mean(w[w>0])
  uWeights<-uWeights/mean(uWeights[uWeights>0])
  weightClusterSize<-as.double(weightClusterSize)
  
  res<-.critFunction(M=M, clu=clu, weights=w, uWeights=uWeights, dimensions=sum(tmNclu), n=n, weightClusterSize=weightClusterSize, diagonal = diagonal, sBorders = limitType, bordersMatLower = bordersMatLower, bordersMatUpper = bordersMatUpper, bordersSeperateLower = bordersSeperateLower, bordersSeperateUpper = bordersSeperateUpper)
  
  # wByHB<-blockmodeling::funByBlocks(w,clu=rep(1:length(n),times=n),ignore.diag=FALSE, FUN=sum)
  # k<-length(unique(clu))
  # parByHB<-findActiveParam(M, n, k, na.rm=TRUE)
  # ICLpen<- sum(unclass(parByHB*log(wByHB)))+sum((k-1)*log(n))
  # ICL<- -res - 1/2*ICLpen

  return(.ICL(M = M,k = k,weights = w,n = n,err=res))

  # res<-list(M=M, clu=clu, IM=IM, err=err, best=list(list(M=M, clu=clu, IM=IM)))
  # return(res)
}

#' @importFrom utils packageVersion

#| Undeclared function for internal use
.ICL<-function(M, k, weights, n, err, ll){
  if(missing(err)) err<- -ll
  w<-weights
  dmW<-dim(w)
  if((packageVersion("blockmodeling")<"1.0.9") & length(dim(w))==3) w<-aperm(w, c(3,1,2))
  wByHB<-blockmodeling::funByBlocks(w,clu=rep(1:length(n),times=n),ignore.diag=FALSE, FUN=sum)
  if((packageVersion("blockmodeling")<"1.0.9") & length(dim(wByHB))==3) wByHB<-aperm(wByHB, c(2,3,1))
  
  if(length(dim(wByHB))==3){
    wByHB<-aperm(wByHB,c(2,3,1))
  }else if(length(dim(wByHB))==2){
    wByHB<-array(wByHB,dim=c(dim(wByHB),1))
  } else { wByHB<-array(wByHB,dim=c(1,1,1))}
  #k<-length(unique(clu))
  parByHB<-findActiveParam(M, n, k, na.rm=TRUE)
  wByHB[wByHB==0]<-1
  ICLpen<- sum(unclass(parByHB*log(wByHB)))+sum((k-1)*log(n))
  ICL<- -err - 1/2*ICLpen
  return(ICL)
}
