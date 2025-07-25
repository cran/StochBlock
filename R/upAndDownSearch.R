#' Perform Up-and-Down Search Optimization
#'
#' The \code{upAndDownSearch} function performs an up-and-down search (increasing and lowering the number of clusters) on a given dataset, starting from an initial partition and fitness value. The function assumes that the higher fitness value is better. Experimental! The \code{stochBlockForUDS} is just a wrapper to \code{\link{stochBlock}} function to be used within \code{upAndDownSearch} function.
#'
#' @importFrom stats runif
#' @import blockmodeling
#'
#' @param data A data object (e.g., data frame or matrix) on which the search is performed.
#' @param initPart Initial partition or configuration to start the search from. It can be a vector for simple datasets and a list of vectors for temporal or linked networks
#' @param initFit Initial fitness value associated with the initial partition.
#' @param optimFun A function used for optimization. A function must accept a data and a initial partition and return a list with an element \code{part} holding the final partition and an element \code{fit} fitness value. 
#' @param nRep Integer. The number of repetitions for the search algorithm. Defaults to 100.
#' @param minPmove Numeric. The minimum probability of moving to a new partition when the fitness does not improve. Defaults to 0.2.
#' @param pFitSumTreshMoveBest Numeric. The threshold for the sum of 1 - probabilities of moving to a new partition before reverting to the best partition. Defaults to 5.
#' @param ... Additional paramters to optimFun.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{data}{The input data.}
#'   \item{finalPart}{The final partition obtained.}
#'   \item{finalFit}{The final fitness value after the search.}
#'   \item{searchHistory}{A list containing the history of partitions and fitness values during the search.}
#'   \item{callUsed}{The call used to invoke the function, capturing the parameters passed.}
#'   \item{initial.param}{A list of initial parameters used in the function call withiut the data.}
#' }
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
#' initClu<-rep(1, times=n)
#' initFit<-ICLStochBlock(M, initClu) # Initial fitness value
#' # Using up-and-down search to optimise the partition
#' res<-upAndDownSearch(data=M,initPart=initClu, initFit=initFit, optimFun=stochBlockForUDS, nRep=10) 
#' plotMat(res$data, clu=res$bestPart) # Have a look at the optimised parition
#' print(res$bestFit) # Print the final fitness value
#' 
#' # Create a synthetic linked-network matrix
#' set.seed(2022)
#' library(blockmodeling)
#' IM<-matrix(c(0.9,.5,0.1,0.8), nrow=2)
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
#' clu1<-rep(1, n)
#' clu2<-rep(2, n)
#' initClu<-list(clu1, clu2)
#' initFit<-ICLStochBlock(Ml, initClu) # Initial fitness value
#' # Using up-and-down search to optimise the partition
#' res<-upAndDownSearch(data=Ml,initPart=initClu, initFit=initFit, optimFun=stochBlockForUDS, nRep=10)
#' plotMat(res$data, clu=res$bestPart) # Have a look at the optimised parition
#' print(res$bestFit) # Print the final fitness value
#' @export
upAndDownSearch<-function(data, initPart, initFit, optimFun, nRep=100, minPmove=0.2, pFitSumTreshMoveBest=5, ...){
  callUsed <- match.call()
  formal_args <- setdiff(names(formals()), "...")
  evaluated_formals <- lapply(formal_args, function(x) eval.parent(substitute(x)))
  names(evaluated_formals) <- formal_args
  
  # Capture evaluated ... arguments
  dot_args <- list(...)
  
  # Combine all initial parameters
  initial.param <- c(evaluated_formals, dot_args)
  initial.param$data<-NULL
  
  searchHistroy<-list()
  linked<-is.list(initPart)
  bestPart<-curPart<-initPart
  bestFit<-curFit<-initFit
  if(linked) {
    sets<-sapply(curPart, length)  
  } else sets<-length(curPart)
  pFitSum<-0
  
  
  on.exit(return(list(data=data, bestPart=bestPart, bestFit=bestFit, searchHistroy=searchHistroy, callUsed=callUsed, initial.param=initial.param)))
  
  for(i in 1:nRep){
    if(i%%2 == 1){ #odd iterations - try increasing number of clusters
      # split one cluster
      if(linked){
        iSet<-sample(length(sets), prob = sets, size = 1)
        while(TRUE){
          newPar<-curPart
          iClu<-sample(newPar[[iSet]],1)
          maxK<-max(newPar[[iSet]])
          iUnits<-which(newPar[[iSet]]==iClu)
          iCluSize<-length(iUnits)
          if(iCluSize>3){
            selUnits<-rbinom(iCluSize, size = 1, prob=0.5)
            if(sum(selUnits)>1 && (iCluSize-sum(selUnits))>1){
              newPar[[iSet]][iUnits[selUnits]]<-maxK+1
              break
            }
          }
        }
      }else{
        while(TRUE){
          newPar<-curPart
          iClu<-sample(newPar,1)
          maxK<-max(newPar)
          iUnits<-which(newPar==iClu)
          iCluSize<-length(iUnits)
          if(iCluSize>3){
            selUnits<-rbinom(iCluSize, size = 1, prob=0.5)==1
            if(sum(selUnits)>1 && (iCluSize-sum(selUnits))>1){
              newPar[iUnits[selUnits]]<-maxK+1
              break
            }
          }
        }
      }
    } else { #even iterations - try decreasing number of clusters
      if(linked){
        for(i2 in 1:3){
          iSet<-sample(length(sets), prob = sets, size = 1)
          if(length(unique(curPart[[iSet]]))<2) next
          newPar<-curPart
          clus<-unique(curPart[[iSet]])
          iClu<-sample(clus,2)
          newPar[[iSet]][newPar[[iSet]]==iClu[2]]<-iClu[1]
          newPar[[iSet]]<-encodeToZeroIndexed(newPar[[iSet]])
          break
        }
      }else{
        if(length(unique(curPart))<2) next
        newPar<-curPart
        clus<-unique(curPart)
        iClu<-sample(clus,2)
        newPar[newPar==iClu[2]]<-iClu[1]
        newPar<-encodeToZeroIndexed(newPar)
      }
    }
    
    
    
    iRes<-optimFun(data, newPar)
    searchHistroy[[i]]<-iRes
    cat(sprintf("Iteration %d/%d\n", i, nRep))
    cat("Number of clusters:\n")
    cat(numClust(iRes$part))
    cat(sprintf("\nNew fitness %f\n", iRes$fit))

    if(iRes$fit>curFit){
      curFit<-iRes$fit
      curPart<-iRes$part
      if(curFit>bestFit){
        cat("Fitness improved!\n")
        bestFit<-curFit
        bestPart<-curPart
        pFitSum<-0
      }
    } else{
      pFit<-exp(iRes$fit-curFit)
      p<-max(pFit, minPmove)
      if(runif(1) <= p){
        curFit<-iRes$fit
        curPart<-iRes$part
        cat("Moved to a worse solution.\n")
      }
      cat("pFit:",pFit,"\n")
      pFitSum<-pFitSum + 1 - pFit
      if(pFitSum > pFitSumTreshMoveBest && pFit!=1){
        pFitSum<-0
        curFit<-bestFit
        curPart<-bestPart
        cat("\nFitness did not improve for a while, reverting to best fitness!\n")
        cat("Number of clusters for currently best solution:\n")
        cat(numClust(bestPart))
        cat("\n")
      }
      
    }
    cat(sprintf("Best fitness %f\n\n", bestFit))
  }
}



#' @rdname upAndDownSearch
#' @export
stochBlockForUDS<-function(data, initPart, ...){
  res<-stochBlock(M=data, clu=initPart, ...)
  return(list(part=res$clu, fit=res$ICL))
}