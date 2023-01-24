#' Computes weights for parts of the multilevel network based on random errors using the SS approach with complete blocks only (compatible with k-means)
#'
#' @import foreach
#' @import parallel
#' @import doParallel
#' @import doRNG
#' @import blockmodeling
#' @importFrom stats na.omit
#' @importFrom stats sd
#' @import Rcpp
#'
#' @param mlNet A multilevel/linked network - The code assumes only one relation --> a matrix.
#' @param cluParts A partition spliting the units into different sets
#' @param k A vecotor of number of clusters for each set of units in the network.
#' @param mWeights The number of repetitions for computing random errors. Defaults to 1000
#' @param sumFun The function to compute the summary of errors, which is then used to compute the weights by computing 1/summary. Defaults to \code{sd}.
#' @param nCores The number of to use for parallel computing. 0 means all available - 1, 1 means only once core - no parallel computing.
#' @param paramGenPar The parameter \code{addParam} from  \code{\link[blockmodeling]{genRandomPar}} (see documentation there). Default here is paramGenPar=list(genPajekPar = FALSE), which is different from the default in \code{\link[blockmodeling]{genRandomPar}}. The same value is used for generating partitions for all partitions.
#' @param weightClusterSize The weight given to cluster sizes. Defalults to 0, as only this is weighted my the tie-based weights.
#' @param ... Paramters passed to \code{\link{llStochBlock}}
#' 
#' @return Weights and "intermediate results":
#'  \item{errArr}{A 3d array of errors (\code{mWeights} for each part of the network)}
#'  \item{errMatSum}{\code{errArr} summed over all repetitions.}
#'  \item{weightsMat}{A matrix of weights, one for each part. An inverse of \code{errMatSum} with NaNs replaced by zeros.}
#'  
#' @references Škulj, D., & Žiberna, A. (2022). Stochastic blockmodeling of linked networks. Social Networks, 70, 240-252.
#' 
#' @author \enc{Aleš, Žiberna}{Ales Ziberna}
#' 
#' @seealso \code{\link{llStochBlock}}; \code{\link{ICLStochBlock}}

weightsMlLoglik<-function(mlNet,cluParts, k, mWeights=1000, sumFun = sd,
                          nCores=0, weightClusterSize=0,
                          paramGenPar=list(genPajekPar = FALSE),...){

  requireNamespace('foreach')
  requireNamespace('doParallel')
  requireNamespace('doRNG')

  if(length(dim(mlNet))!=2){
    stop("Currently the function only work if mlNet is (2-dimensional) matrix!")
    }
  if(nCores == 0){
    nCores <- detectCores()-1
  }

  if(nCores>1 & !getDoParRegistered()){
    registerDoParallel(nCores)
  }

  # if(requireNamespace('blockmodeling')){
    pack<-"blockmodeling"
#  }else{
    # requireNamespace("blockmodelingTest")
    # pack<-"blockmodelingTest"
#  }
  pack<-c(pack,"StochBlock")
  cluParts<-as.numeric(factor(cluParts))
  nn<-table(cluParts)
  mlOrNet<-apply(mlNet, c(1,2),sum)
  parts<-fun.by.blocks(mlOrNet, clu = cluParts, ignore.diag = FALSE, FUN = ss)>0
  errArr<-array(NA,dim=c(dim(parts),mWeights))
  errUnitsMat<-matrix(NA,ncol=length(k),nrow=mWeights)
  
  for(i1 in 1:dim(parts)[1]){
    for(i2 in 1:dim(parts)[2]){
      if(is.finite(parts[i1,i2])&parts[i1,i2]){
        net<-mlOrNet[cluParts==i1, cluParts==i2]
        if(i1!=i2){
          tnet<-matrix(0, nrow=(nn[i1] + nn[i2]), ncol=(nn[i1] + nn[i2]))
          tnet[cluParts[cluParts %in% c(i1,i2)]==i1, cluParts[cluParts %in% c(i1,i2)]==i2]<-net
          net<-tnet
        }
        tmp<-foreach(i =1:mWeights,.combine = rbind, .packages =  pack)%dorng%{
          if(i1==i2){
            tclu<-genRandomPar(k = k[i1],n = nn[i1],addParam=paramGenPar)
            tabClu<-table(tclu)
            ptabClu<-tabClu/nn[i1]
            lptabClu<-log(ptabClu)
            errUnits<- -sum(tabClu*lptabClu)
          } else {
            errUnits<-NA
            tclu<-genRandomPar(k = k[c(i1,i2)],n = nn[c(i1,i2)],addParam=paramGenPar)
# the bellow lines are not needed due to improvment in llStochBlock
#            tclu[[2]]<-tclu[[2]]+k[i1]
#            tclu<-unlist(tclu)
          }
          tCF<-llStochBlock(net, clu=tclu, weightClusterSize=weightClusterSize, ...)
          c(tCF,errUnits)
        }
        errArr[i1,i2,]<-tmp[,1]
        if(i1==i2) errUnitsMat[,i1]<-tmp[,2]
      }
      cat("i1 = ", i1, ", i2 = ", i2,"\n")
    }
  }
  errMatSum<-apply(errArr,c(1,2),sumFun)
  weightsMat<-relInv2(errMatSum)
  errUnitsSum<-apply(errUnitsMat,2,sumFun)
  weightsUnits<-relInv2(errUnitsSum)
  res<-list(ties=list(errArr=errArr, errMatSum=errMatSum, weightsMat=weightsMat),units=list(errMat=errUnitsMat, errSum=errUnitsSum, weights=weightsUnits))
  return(res)
}

