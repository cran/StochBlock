#| Un-exported function for internal use from blockmodeling
unlistPar<-function(part){
  if (is.list(part)) {
    part <- sapply(part, paste, collapse = " ")
    part <- paste(paste("\nMode ", 1:length(part), ":", sep = ""),
                  part, collapse = "", sep = " ")
  }
  part
}

#' Function that performs stochastic one-mode and linked blockmodeling by optimizing a single partition. If \code{clu} is a list, the method for linked/multilevel networks is applied
#'
#' @import doParallel
#' @import doRNG
#' @import blockmodeling
#'
#' @param M A matrix representing the (usually valued) network. For multi-relational networks, this should be an array with the third dimension representing the relation.
#' @param clu A partition. Each unique value represents one cluster. If the network is one-mode, than this should be a vector, else a list of vectors, one for each mode. Similarly, if units are comprised of several sets, clu should be the list containing one vector for each set.
#' @param weights The weights for each cell in the matrix/array. A matrix or an array with the same dimensions as \code{M}.
#' @param uWeights The weights for each unin. A vector with the length equal to the number of units (in all sets).
#' @param diagonal How should the diagonal values be treated. Possible values are:
#' \itemize{
#'   \item ignore - diagonal values are ignored
#'   \item seperate - diagonal values are treated seperately
#'   \item same - diagonal values are treated the same as all other values
#' }
#' @param limits If \code{diagonal} is \code{"ignore"} or \code{"same"}, an array with dimensions equal to:
#' \itemize{
#'   \item number of clusters (of all types)
#'   \item number of clusters (of all types)
#'   \item number of relations
#'   \item 2 - the first is lower limit and the second is upper limit
#' }
#' If \code{diagonal} is \code{"seperate"}, a list of two array. The first should be as described above, representing limits for off diagonal values. The second should be similar with only 3 dimensions, as one of the first two must be omitted.
#' @param limitType Type of limit to use. Forced to 'none' if \code{limits} is \code{NULL}. Otherwise, one of either \code{outer} or \code{inner}.
#' @param weightClusterSize The weight given to cluster sizes (logprobabilites) compared to ties in loglikelihood. Defaults to 1, which is "classical" stochastic blockmodeling.
#' @param addOne Should one tie with the value of the tie equal to the density of the superBlock be added to each block to prevent block means equal to 0 or 1 and also "shrink" the block means toward the superBlock mean. Defaults to TRUE.
#' @param eps If addOne = FALSE, the minimal deviation from 0 or 1 that the block mean/density can take.
#'
#' @return A list of class \code{opt.par} normally passed other commands with \code{StockBlockORP} and containing:
#'\item{clu}{A vector (a list for multi-mode networks) indicating the cluster to which each unit belongs;}
#'\item{IM}{Image matrix of this partition;}
#'\item{weights}{The weights for each cell in the matrix/array. A matrix or an array with the same dimensions as \code{M}.}
#'\item{uWeights}{The weights for each unit. A vector with the length equal to the number of units (in all sets).}
#'\item{err}{The error as the sum of the inconsistencies between this network and the ideal partitions.}
#'\item{ICL}{Integrated Criterion Likelihood for this partition}
#'
#' @references Škulj, D., & Žiberna, A. (2022). Stochastic blockmodeling of linked networks. Social Networks, 70, 240-252.
#' 
#' @author \enc{Aleš, Žiberna}{Ales Ziberna}
#' 
#' @seealso \code{\link{stochBlockORP}}
#' 
#' @examples
#' 
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
#' res<-stochBlock(M,clu) # Optimising the partition
#' plot(res) # Have a look at the optimised parition
#' 
#' # Create a synthetic linked-network matrix
#' set.seed(2022)
#' library(blockmodeling)
#' IM<-matrix(c(0.8,.4,0.2,0.8), nrow=2)
#' clu<-rep(1:2, each=20) # Partition to generate
#' n<-length(clu)
#' nClu<-length(unique(clu)) # Number of clusters to generate
#' M1<-matrix(rbinom(n^2,1,IM[clu,clu]),ncol=n, nrow=n) # First network
#' M2<-matrix(rbinom(n^2,1,IM[clu,clu]),ncol=n, nrow=n) # Second network
#' M12<-diag(n) # Linking network
#' nn<-c(n,n)
#' k<-c(2,2)
#' Ml<-matrix(0, nrow=sum(nn),ncol=sum(nn)) 
#' Ml[1:n,1:n]<-M1
#' Ml[n+1:n,n+1:n]<-M2
#' Ml[n+1:n, 1:n]<-M12 
#' plotMat(Ml) # Linked network
#' clu1<-sample(1:2,nrow(M1),replace=TRUE)
#' clu2<-sample(3:4,nrow(M1),replace=TRUE)
#' plotMat(Ml,list(clu1,clu2)) # Have a look at this random partition
#' res<-stochBlock(Ml,list(clu1,clu2)) # Optimising the partition
#' plot(res) # Have a look at the optimised parition
#' 
#' @useDynLib StochBlock, .registration = TRUE
#' @export

stochBlock<-function(M,
                  clu,
                  weights=NULL,
				  uWeights=NULL,
				  diagonal = c("ignore","seperate","same"),
				  limitType=c("none","inside","outside"),
                  limits=NULL,
				  weightClusterSize = 1,
				  addOne = TRUE,
				  eps = 0.001){
  
  if(!all(unique(as.vector(unlist(unclass(M))))%in%c(0,1))) stop("Data must be binary (only 0 and 1)!")
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
  if(is.null(uWeights)){
	uWeights<-rep(1.0, n1)
  }
  if(length(uWeights)!=n1) stop("uWeights has wrong length!")
  w<-weights

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
  	  bordersSeperateLower <- bordersSeperateUpper<-NULL
  	  if(is.list(limits)){
  	    limits<-limits[[1]]
  	  }
  	  if(!is.array(limits)) stop("'limits' must be specified as an array!")
  	  dl<-dim(limits)
  	  if(length(dl)==3 & dim(M)[3]==1) limits<-array(limits, dim=c(dl[1:2],1,dl[3]))
	  dl<-dim(limits)
  	  if(length(dl)!=4) stop("'limits' has wrong dimensions (see help for correct dimensions)")
  	  
  	  if(all(dim(limits)!=c(sum(tmNclu),sum(tmNclu),dim(M)[3],2))){
  	    stop("'limits' has wrong dimensions (see help for correct dimensions)")
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
  	  dl<-dim(limits)
	  if(length(dl)!=4) stop("First element of 'limits' has wrong dimensions (see help for correct dimensions)")
  	  if(all(dim(limits)!=c(sum(tmNclu),sum(tmNclu),dim(M)[3],2))){
  	    stop("First element of 'limits' has wrong dimensions (see help for correct dimensions)")
  	  } else{
  	    bordersMatLower <- limits[,,,1]
  	    bordersMatUpper <- limits[,,,2]
  	  }
  
  	  if(!is.array(diagLimits)) stop("Second element of 'limits' must be specified as an array!")
  	  dl<-dim(diagLimits)
  	  if(length(dl)==2 & dim(M)[3]==1) limits<-array(limits, dim=c(dl[1],1,dl[2]))
	  dl<-dim(diagLimits)
  	  if(length(dl)!=3) stop("Second element of 'limits' has wrong dimensions (see help for correct dimensions)")
  	  if(all(dim(diagLimits)!=c(sum(tmNclu),dim(M)[3],2))){
  	    stop("Second element of 'limits' has wrong dimensions (see help for correct dimensions)")
  	  } else{
  	    bordersSeperateLower <- diagLimits[,,1]
  	    bordersSeperateUpper <- diagLimits[,,2]
  	  }
  	}
  }
  
  if(diagonal == "ignore")for(i in 1:dim(w)[3]){
    diag(w[,,i])<-0
  } 

  #normalization of weights
  w<-w*findEmptySuperbocks(M,n = n)
  w<-w/mean(w[w>0])
  uWeights<-uWeights/mean(uWeights[uWeights>0])


  weightClusterSize<-as.double(weightClusterSize)
  
  res<-.kmBlock(M=M, clu=clu, weights=w, uWeights=uWeights, n=n, nClu=tmNclu, diagonal = diagonal, weightClusterSize = weightClusterSize,  sBorders = limitType, bordersMatLower = bordersMatLower, bordersMatUpper = bordersMatUpper, bordersSeperateLower = bordersSeperateLower, bordersSeperateUpper = bordersSeperateUpper, addOne = addOne, eps = eps)
  
	  
  res<-list(M=M, clu=blockmodeling::splitClu(res$bestClu,n), IM=res$IM, err=res$bestCf, weights=w, uWeights=uWeights, n=n, ICL=.ICL(M=M, k = k, weights=w, n=n, err=res$bestCf))
  #return(res)
  class(res)<-"opt.par"
  return(res)
}





#' Function that computes criterion function used in stochastic one-mode and linked blockmodeling. If \code{clu} is a list, the method for linked/multilevel networks is applied
#'
#' @param M A matrix representing the (usually valued) network. For multi-relational networks, this should be an array with the third dimension representing the relation.
#' @param clu A partition. Each unique value represents one cluster. If the network is one-mode, than this should be a vector, else a list of vectors, one for each mode. Similarly, if units are comprised of several sets, clu should be the list containing one vector for each set.
#' @param weights The weights for each cell in the matrix/array. A matrix or an array with the same dimensions as \code{M}.
#' @param uWeights The weights for each unit. A vector with the length equal to the number of units (in all sets).
#' @param diagonal How should the diagonal values be treated. Possible values are:
#' \itemize{
#'   \item ignore - diagonal values are ignored
#'   \item seperate - diagonal values are treated separately
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
#' @param weightClusterSize The weight given to cluster sizes (log-probabilities) compared to ties in loglikelihood. Defaults to 1, which is "classical" stochastic blockmodeling.
#'
#' @return - the value of the log-likelihood criterion for the partition \code{clu} on the network represented by \code{M} for binary stochastic blockmodel.
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
#' ll_pre<-llStochBlock(M,clu) # Calculate its loglikelihood
#' res<-stochBlockORP(M,k=2,rep=10) # Optimizing the partition
#' plot(res) # Have a look at the optimized partition
#' ll_post<-llStochBlock(M,clu(res)) # Calculate its loglikelihood
#' # We expect the loglikelihood pre-optimization to be smaller:
#' (-ll_pre)<(-ll_post)
#'   
#'
#' @export

llStochBlock<-function(M,
                       clu,
                       weights=NULL,
                       uWeights=NULL,
                       diagonal = c("ignore","seperate","same"),
                       limitType=c("none","inside","outside"),
                       limits=NULL,
                       weightClusterSize=1.0,
                       addOne = TRUE,
                       eps = 0.001){

  if(!all(unique(as.vector(unlist(unclass(M))))%in%c(0,1))) stop("Data must be binary (only 0 and 1)!")
  n1<-dim(M)[1]
  if(is.list(clu)) {
    n<-sapply(clu, length)
  }  else{
    n<-length(clu)
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
      if(dim(limits)!=4) stop("'limits' has wrong dimensions (see help for correct dimensions)")
      
      if(all(dim(limits)!=c(sum(tmNclu),sum(tmNclu),dim(M)[3],2))){
        stop("'limits' has wrong dimensions (see help for correct dimensions)")
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
      if(dim(limits)!=4) stop("First element of 'limits' has wrong dimensions (see help for correct dimensions)")
      if(all(dim(limits)!=c(sum(tmNclu),sum(tmNclu),dim(M)[3],2))){
        stop("First element of 'limits' has wrong dimensions (see help for correct dimensions)")
      } else{
        bordersMatLower <- limits[,,,1]
        bordersMatUpper <- limits[,,,2]
      }
      
      if(!is.array(diagLimits)) stop("Second element of 'limits' must be specified as an array!")
      dl<-dim(diagLimits)
      if(length(dl)==2 & dim(M)[3]==1) limits<-array(limits, dim=c(dl[1],1,dl[2]))
      if(dim(diagLimits)!=3) stop("Second element of 'limits' has wrong dimensions (see help for correct dimensions)")
      if(all(dim(diagLimits)!=c(sum(tmNclu),dim(M)[3],2))){
        stop("Second element of 'limits' has wrong dimensions (see help for correct dimensions)")
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
  
  res<-.critFunction(M=M, clu=clu, weights=w, uWeights=uWeights, dimensions=sum(tmNclu), n=n, weightClusterSize=weightClusterSize, diagonal = diagonal, sBorders = limitType, bordersMatLower = bordersMatLower, bordersMatUpper = bordersMatUpper, bordersSeperateLower = bordersSeperateLower, bordersSeperateUpper = bordersSeperateUpper, addOne = addOne, eps = eps)
  return(res)
  
  # res<-list(M=M, clu=clu, IM=IM, err=err, best=list(list(M=M, clu=clu, IM=IM)))
  # return(res)
}

#' A function for optimizing multiple random partitions using stochastic one-mode and linked blockmodeling. Similar to optRandomParC, but calling stochBlock for optimizing individual partitions.
#'
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import doRNG
#' @import blockmodeling
#' @importFrom stats var
#' @importFrom stats rbinom
#' @importFrom stats kmeans
#'
#' @param M A square matrix giving the adjaciency relationg between the network's nodes (aka vertexes)
#' @param k The number of clusters used in the generation of partitions.
#' @param rep The number of repetitions/different starting partitions to check.
#' @param save.initial.param Should the inital parameters(\code{approaches}, ...) of using \code{stochBlock} be saved. The default value is \code{TRUE}.
#' @param deleteMs Delete networks/matrices from the results of to save space. Defaults to \code{TRUE}.
#' @param max.iden Maximum number of results that should be saved (in case there are more than \code{max.iden} results with minimal error, only the first \code{max.iden} will be saved).
#' @param return.all If \code{FALSE}, solution for only the best (one or more) partition/s is/are returned.
#' @param return.err Should the error for each optimized partition be returned. Defaults to \code{TRUE}.
#' @param seed Optional. The seed for random generation of partitions.
#' @param parGenFun The function (object) that will generate random partitions. The default function is   \code{\link{genRandomPar}}. The function has to accept the following parameters: \code{k} (number o of partitions by modes, \code{n} (number of units by modes), \code{seed} (seed value for random generation of partition), \code{addParam} (a list of additional parameters).
#' @param mingr Minimal allowed group size.
#' @param maxgr Maximal allowed group size.
#' @param addParam A list of additional parameters for function specified above. In the usage section they are specified for the default function \code{\link{genRandomPar}}.
#' @param maxTriesToFindNewPar The maximum number of partition try when trying to find a new partition to optimize that was not yet checked before - the default value is \code{rep * 1000}.
#' @param skip.par The partitions that are not allowed or were already checked and should therefore be skipped.
#' @param printRep Should some information about each optimization be printed.
#' @param n The number of units by "modes". It is used only for generating random partitions. It has to be set only if there are more than two modes or if there are two modes, but the matrix representing the network is one mode (both modes are in rows and columns).
#' @param nCores Number of cores to be used. Value \code{0} means all available cores. It can also be a cluster object.
#' @param useParLapply Should \code{parLapplyLB} be used (otherwise \code{foreach} is used). Defaults to true as it needs less dependencies. It might be removed in future releases and only allow the use of parLapplyLB.
#' @param cl The cluster to use (if formed beforehand). Defaults to \code{NULL}.
#' @param stopcl Should the cluster be stopped after the function finishes. Defaults to \code{is.null(cl)}.
#' @param \dots Arguments passed to other functions, see \code{\link{stochBlock}}.
#'
#' @return A list of class "opt.more.par" containing:
#'  \item{M}{The one- or multi-mode matrix of the network analyzed}
#'   \item{res}{If \code{return.all = TRUE} - A list of results the same as \code{best} - one \code{best} for each partition optimized.}
#'   \item{best}{A list of results from \code{stochblock}, only without \code{M}.}
#'   \item{err}{If \code{return.err = TRUE} - The vector of errors or inconsistencies = -log-likelihoods.}
#'   \item{ICL}{Integrated classification likelihood for the best partition.}
#'   \item{checked.par}{If selected - A list of checked partitions. If \code{merge.save.skip.par} is \code{TRUE}, this list also includes the partitions in \code{skip.par}.}
#'   \item{call}{The call to this function.}
#'   \item{initial.param}{If selected - The initial parameters are used.}
#'   \item{Random.seed}{.Random.seed at the end of the function.}
#'   \item{cl}{Cluster used for parallel computations if supplied as an input parameter.}
#'   
#' @section Warning:
#' It should be noted that the time needed to optimise the partition depends on the number of units (aka nodes) in the networks as well as the number of clusters
#' due to the underlying algorithm. Hence, partitioning networks with 100 units and large number of blocks (e.g., >5) can take a long time (from 20 minutes to a few hours or even days).
#' 
#' @references Škulj, D., & Žiberna, A. (2022). Stochastic blockmodeling of linked networks. Social Networks, 70, 240-252.
#' 
#' @author \enc{Aleš, Žiberna}{Ales Ziberna}
#'
#' @examples
#'# Simple one-mode network
#'library(blockmodeling)
#'k<-2
#'blockSizes<-rep(20,k)
#'IM<-matrix(c(0.8,.4,0.2,0.8), nrow=2)
#'if(any(dim(IM)!=c(k,k))) stop("invalid dimensions")
#'
#'set.seed(2021)
#'clu<-rep(1:k, times=blockSizes)
#'n<-length(clu)
#'M<-matrix(rbinom(n*n,1,IM[clu,clu]),ncol=n, nrow=n)
#'diag(M)<-0
#'plotMat(M)
#'
#'resORP<-stochBlockORP(M,k=2, rep=10, return.all = TRUE)
#'resORP$ICL
#'plot(resORP)
#'clu(resORP)
#'
#'
#'# Linked network
#'library(blockmodeling)
#'set.seed(2021)
#'IM<-matrix(c(0.8,.4,0.2,0.8), nrow=2)
#'clu<-rep(1:2, each=20)
#'n<-length(clu)
#'nClu<-length(unique(clu))
#'M1<-matrix(rbinom(n^2,1,IM[clu,clu]),ncol=n, nrow=n)
#'M2<-matrix(rbinom(n^2,1,IM[clu,clu]),ncol=n, nrow=n)
#'M12<-diag(n)
#'nn<-c(n,n)
#'k<-c(2,2)
#'Ml<-matrix(0, nrow=sum(nn),ncol=sum(nn))
#'Ml[1:n,1:n]<-M1
#'Ml[n+1:n,n+1:n]<-M2
#'Ml[n+1:n, 1:n]<-M12
#'plotMat(Ml)
#'
#'resMl<-stochBlockORP(M=Ml, k=k, n=nn, rep=10)
#'resMl$ICL
#'plot(resMl)
#'clu(resMl)
#'
#' @author \enc{Aleš, Žiberna}{Ales Ziberna}
#' 
#' @export

stochBlockORP<-function(M, #a square matrix
                         k,#number of clusters/groups
                         rep,#number of repetitions/different starting partitions to check
                         save.initial.param=TRUE,  #save the initial parameters of this call
                         deleteMs=TRUE, #delete networks/matrices from results of optParC or optParMultiC to save space
                         max.iden=10, #the maximum number of results that should be saved (in case there are more than max.iden results with minimal error, only the first max.iden will be saved)
                         return.all=FALSE,#if 'FALSE', solution for only the best (one or more) partition/s is/are returned
                         return.err=TRUE,#if 'FALSE', only the results of crit.fun are returned (a list of all (best) solutions including errors), else the result is list
                         seed=NULL,#the seed for random generation of partitions
                         parGenFun = blockmodeling::genRandomPar, #The function that will generate random partitions. It should accept arguments: k (number of partitions by modes, n (number of units by modes), seed (seed value for random generation of partition), addParam (a list of additional parameters)
                         mingr=NULL, #minimal allowed group size (defaults to c(minUnitsRowCluster,minUnitsColCluster) if set, else to 1) - only used for parGenFun function 
                         maxgr=NULL, #maximal allowed group size (default to c(maxUnitsRowCluster,maxUnitsColCluster) if set, else to Inf) - only used for parGenFun function 
                         addParam=list(  #list of additional parameters for generating partitions. Here they are specified for the default function "genRandomPar"
                           genPajekPar = TRUE,     #Should the partitions be generated as in Pajek (the other options is completely random)
                           probGenMech = NULL),    #Here the probabilities for different mechanisms for specifying the partitions are set. If not set this is determined based on the previous parameter.
                         maxTriesToFindNewPar=rep*10,    #The maximum number of partition try when trying to find a new partition to optimize that was not yet checked before 
                         skip.par = NULL, #partitions to be skipped
                         printRep= ifelse(rep<=10,1,round(rep/10)), #should some information about each optimization be printed
                         n=NULL, #the number of units by "modes". It is used only for generating random partitions. It has to be set only if there are more than two modes or if there are two modes, but the matrix representing the network is onemode (both modes are in rows and columns)
                         nCores=1, #number of cores to be used 0 -means all available cores, can also be a cluster object,
                         useParLapply=FALSE, #should parLapply be used instead of foreach
                         cl = NULL, #the cluster to use (if formed beforehand)
                         stopcl = is.null(cl), # should the cluster be stopped
                         ... #paramters to stochBlock
 ){
   dots<-list(...)
   
   if(!all(unique(as.vector(unlist(unclass(M))))%in%c(0,1))) stop("Data must be binary (only 0 and 1)!")
   
   if(save.initial.param)initial.param<-c(tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")),dots=list(...))#saves the inital parameters
   
   if(is.null(mingr)){
     if(is.null(dots$minUnitsRowCluster)){
       mingr<-1
     } else {
       mingr<-c(dots$minUnitsRowCluster,dots$minUnitsColCluster)
     }
   }
   
   if(is.null(maxgr)){
     if(is.null(dots$maxUnitsRowCluster)){
       maxgr<-Inf
     } else {
       maxgr<-c(dots$maxUnitsRowCluster,dots$maxUnitsColCluster)
     }
   }
   
   nmode<-length(k)
   
   res<-list(NULL)
   err<-NULL
   dots<-list(...)
   
   if(save.initial.param)initial.param<-c(tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")),dots=list(...))#saves the inital parameters
   
   
   if(is.null(n)) if(nmode==1){
     n<-dim(M)[1]
   } else if(nmode==2){
     n<-dim(M)[1:2]
   } else warning("Number of nodes by modes can not be determined. Parameter 'n' must be supplied!!!")
   
   if(!is.null(seed))set.seed(seed)
   
   
   on.exit({
     res1 <- res[which(err==min(err, na.rm = TRUE))]
     best<-NULL
     best.clu<-NULL
     for(i in 1:length(res1)){
       for(j in 1:length(res1[[i]]$best)){
         if(
           ifelse(is.null(best.clu),
                  TRUE,
                  if(nmode==1){
                    !any(sapply(best.clu,rand2,clu2=res1[[i]]$clu)==1)
                  } else {
                    !any(sapply(best.clu,function(x,clu2)rand2(unlist(x),clu2),clu2=unlist(res1[[i]]$clu))==1)
                   }
           )
         ){
           best<-c(best,res1[i])
           best.clu<-c(best.clu,list(res1[[i]]$clu))
         }
         
         if(length(best)>=max.iden) {
           warning("Only the first ",max.iden," solutions out of ",length(na.omit(err))," solutions with minimal -loglikelihood will be saved.\n")
           break
         }
         
       }
     }
     
     names(best)<-paste("best",1:length(best),sep="")
     
     if(any(na.omit(err)==-Inf) || ss(na.omit(err))!=0 || length(na.omit(err))==1){
       cat("\n\nOptimization of all partitions completed\n")
       cat(length(best),"solution(s) with minimal -loglikelihood =", min(err,na.rm=TRUE), "found.","\n")
     }else {
       cat("\n\nOptimization of all partitions completed\n")
       cat("All",length(na.omit(err)),"solutions have -loglikelihood",err[1],"\n")
     }
     
     call<-list(call=match.call())
     best<-list(best=best)
     checked.par<-list(checked.par=skip.par)
     if(return.all) res<-list(res=res) else res<-NULL
     if(return.err) err<-list(err=err) else err<-NULL
     if(!exists("initial.param")){
       initial.param<-NULL
     } else initial.param=list(initial.param)
     
     res<-c(list(M=M),list(ICL=best[[1]][[1]]$ICL),res,best,err,checked.par,call,initial.param=initial.param, list(Random.seed=.Random.seed, cl=cl))
     class(res)<-"opt.more.par"
     return(res)
   })
   
   
   
   if(nCores==1||!requireNamespace('parallel')){
     if(nCores!=1) {
       warning("Only single core is used as package 'parallel' is not available", immediate.=TRUE)
     }
     for(i in 1:rep){
       if(printRep & (i%%printRep==0)) cat("\n\nStarting optimization of the partiton",i,"of",rep,"partitions.\n")
       find.unique.par<-TRUE
       ununiqueParTested=0
       while(find.unique.par){
         temppar<-parGenFun(n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam)
         
         find.unique.par<-
           ifelse(is.null(skip.par),
                  FALSE,
                  if(nmode==1) {
                    any(sapply(skip.par,rand2,clu2=temppar)==1)
                  } else any(sapply(skip.par,function(x,clu2)rand2(unlist(x),clu2),clu2=unlist(temppar))==1)
           )
         ununiqueParTested=ununiqueParTested+1
         endFun<-ununiqueParTested>=maxTriesToFindNewPar
         if(endFun) {
           break
         } else if(ununiqueParTested%%10==0) cat(ununiqueParTested,"partitions tested for unique partition\n")
       }
       
       if(endFun) break
       
       skip.par<-c(skip.par,list(temppar))
       
       if(printRep==1) cat("Starting partition:",unlistPar(temppar),"\n")
       res[[i]]<-stochBlock(M=M, clu=temppar,  ...)
       if(deleteMs){
         res[[i]]$M<-NULL
       }
       res[[i]]$best<-NULL
 
       err[i]<-res[[i]]$err
       if(printRep==1) cat("Final -loglikelihood:",err[i],"\n")
       if(printRep==1) cat("Final partition:   ",unlistPar(res[[i]]$clu),"\n")
     }
   } else {
     oneRep<-function(i,M,n,k,mingr,maxgr,addParam,rep, parGenFun,...){
       temppar<-parGenFun(n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam)
       #skip.par<-c(skip.par,list(temppar))
       
       tres <- try(stochBlock(M=M, clu=temppar,  ...))
       if(inherits(x = tres,what = "try-error")){
         tres<-list("try-error"=tres, err=Inf, startPart=temppar)
       }
       if(deleteMs){
         tres$M<-NULL
       }
       tres$best<-NULL
       return(list(tres))
     }
     
     if(!requireNamespace("doParallel")|!requireNamespace("doRNG")) useParLapply<-TRUE
     
     if(nCores==0){
       nCores<-detectCores()-1                    
     }
     
 	pkgName<-utils::packageName()
 	if(is.null(pkgName)) pkgName<-utils::packageName(environment(fun.by.blocks))
     if(useParLapply) {
       if(is.null(cl)) cl<-makeCluster(nCores)
       clusterSetRNGStream(cl)
       nC<-nCores
       #clusterExport(cl, varlist = c("kmBlock","kmBlockORP"))
       #clusterExport(cl, varlist = "kmBlock")
       clusterExport(cl, varlist = "pkgName", envir=environment())
       clusterEvalQ(cl, expr={requireNamespace(pkgName,character.only = TRUE)})
       res<-parLapplyLB(cl = cl,1:rep, fun = oneRep, M=M,n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam,rep=rep, parGenFun=parGenFun,...)
       if(stopcl) stopCluster(cl)
       res<-lapply(res,function(x)x[[1]])
     } else {
       requireNamespace("doParallel")
       requireNamespace("doRNG")
       if(!getDoParRegistered()|(getDoParWorkers()!=nCores)){
         if(!is.null(cl)) {
           #cl<-makeCluster(nCores)
           registerDoParallel(cl)
         } else registerDoParallel(nCores)
       }
       nC<-getDoParWorkers()
       
       res<-foreach(i=1:rep,.combine=c, .packages=pkgName) %dorng% oneRep(i=i,M=M,n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam,rep=rep, parGenFun=parGenFun,...)
       if(!is.null(cl) & stopcl) {
         registerDoSEQ()
         stopCluster(cl)
       }
     }
 	err<-sapply(res,function(x)x$err)    
   }
}


# Internal use
findEmptySuperbocks<-function(M, n, na.rm=TRUE){
  if(length(n)==1) return(1)
  if(sum(n)!=dim(M)[1]) stop("Dimensions do not match!")
  kVec<-1:length(n)
  clu<-rep(kVec,times=n)
  w<-M
  w[]<-1
  for(r in 1:dim(M)[3]){
    for(i in kVec)for(j in kVec){
      w[clu==i, clu==j,r]<-var(as.vector(M[clu==i, clu==j,r]), na.rm = na.rm)>0
    }  
  }
  return(w)
}

