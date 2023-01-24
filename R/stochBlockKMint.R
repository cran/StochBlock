#' A function for using k-means to initialized the  stochastic one-mode and linked blockmodeling.
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
#' @inheritParams stochBlockORP
#' @param nstart number of random starting points for the classical k-means algorithm (for each set of units). Defaults to \code{100}.
# @param iter.max Maximum number of iterations from each random starting point (according to \code{nstart}). Defaults to \code{1000}.
#' @param perm Number or partitions obtained by randomly permuting the k-means partition - if 0, no permutations are made, only the original partition is analyzed.
#' @param sharePerm The probability that a unit will have their randomly assigned. Defaults to \code{0.20}.
#'
#' @return A list containing:
#'  \item{M}{The one- or multi-mode matrix of the network analyzed}
#'   \item{res}{If \code{return.all = TRUE} - A list of results the same as \code{best} - one \code{best} for each partition optimized.}
#'   \item{best}{A list of results from \code{stochblock}, only without \code{M}.}
#'   \item{err}{If \code{return.err = TRUE} - The vector of errors or inconsistencies of the empirical  network with the ideal partitions.}
#'   \item{nIter}{The vector of the iterations on each starting partition. If many of the values equal\code{maxiter}, then  \code{maxiter} may be too small.}
#'   \item{checked.par}{If selected - A list of checked partitions. If \code{merge.save.skip.par} is \code{TRUE}, this list also includes the partitions in \code{skip.par}.}
#'   \item{call}{The call to this function.}
#'   \item{initial.param}{If selected - The initial parameters are used.}
#'   
#' 
#' @references Škulj, D., & Žiberna, A. (2022). Stochastic blockmodeling of linked networks. Social Networks, 70, 240-252.
#' 
#' @author \enc{Aleš, Žiberna}{Ales Ziberna}
#'
#' @export

stochBlockKMint<-function(M, #a square matrix
                        k,#number of clusters/groups
                        nstart = 100, #number of random starting points for the classical k-means algorithm (for each set of units).
                        perm = 0, #number or partitions obtained by randomly permuting the k-means partition - if 0, no permutations are made, only the original partition is analyzed. 
                        sharePerm = 0.20, #the probability that a unit will have their randomly assigned. 
                        save.initial.param=TRUE,  #save the initial parametrs of this call
                        deleteMs=TRUE, #delete networks/matrices from results of optParC or optParMultiC to save space
                        max.iden=10, #the maximum number of results that should be saved (in case there are more than max.iden results with minimal error, only the first max.iden will be saved)
                        return.all=FALSE,#if 'FALSE', solution for only the best (one or more) partition/s is/are returned
                        return.err=TRUE,#if 'FALSE', only the resoults of crit.fun are returned (a list of all (best) soulutions including errors), else the resoult is list
                        seed=NULL,#the seed for random generation of partitions
#                        mingr=NULL, #minimal allowed group size (defaults to c(minUnitsRowCluster,minUnitsColCluster) if set, else to 1) - only used for parGenFun function 
#                        maxgr=NULL, #maximal allowed group size (default to c(maxUnitsRowCluster,maxUnitsColCluster) if set, else to Inf) - only used for parGenFun function 
                        maxTriesToFindNewPar=perm*10,    #The maximum number of partition try when trying to find a new partition to optimize that was not yet checked before 
                        skip.par = NULL, #partitions to be skiped
                        printRep= ifelse(perm<=10,1,round(perm/10)), #should some information about each optimization be printed
                        n=NULL, #the number of units by "modes". It is used only for generating initial partitions. It has to be set only if there are more than two modes or if there are two modes, but the matrix representing the network is onemode (both modes are in rows and columns)
                        nCores=1, #number of cores to be used 0 -means all available cores, can also be a cluster object,
                        useParLapply=FALSE, #should parLapply be used instead of foreach
                        cl = NULL, #the cluster to use (if formed beforehand)
                        stopcl = is.null(cl), # should the cluster be stoped
                        ... #paramters to stochBlock
){
  dots<-list(...)
  
  if(save.initial.param)initial.param<-c(tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")),dots=list(...))#saves the inital parameters
  if(nCores!=1) {
    warning("Mullticore implementation is not yet supported!")
    nCores<-1
  }
  # if(is.null(mingr)){
  #   if(is.null(dots$minUnitsRowCluster)){
  #     mingr<-1
  #   } else {
  #     mingr<-c(dots$minUnitsRowCluster,dots$minUnitsColCluster)
  #   }
  # }
  # 
  # if(is.null(maxgr)){
  #   if(is.null(dots$maxUnitsRowCluster)){
  #     maxgr<-Inf
  #   } else {
  #     maxgr<-c(dots$maxUnitsRowCluster,dots$maxUnitsColCluster)
  #   }
  # }
  
  nmode<-length(k)
  
  res<-list(NULL)
  err<-NULL
  dots<-list(...)
  
  if(save.initial.param)initial.param<-c(tryCatch(lapply(as.list(sys.frame(sys.nframe())),eval),error=function(...)return("error")),dots=list(...))#saves the inital parameters
  
  
  if(is.null(n)) if(nmode==1){
    n<-dim(M)[1]
  } else if(nmode==2){
    n<-dim(M)[1:2]
  } else stop("Number of nodes by modes can not be determined. Parameter 'n' must be supplied!!!")
  
  if(nmode!=length(n)) stop("The lengths of k and n must match!")
    
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
  
  dat<-cbind(M, t(M))
  modeVec<-rep(1:nmode, n)
  part<-list()

  for(iMode in 1:nmode){
    iDat<-dat[modeVec==iMode,]
    iDat<-iDat[,colSums(iDat^2)>0]
    iKmRes<-kmeans(iDat, centers = k[iMode],nstart = nstart)
    part[[iMode]]<-iKmRes$cluster
  }
  perm <- perm+1
  
  if(nCores==1||!requireNamespace("parallel")){
    if(nCores!=1) {
      warning("Only single core is used as package 'parallel' is not available", immediate.=TRUE)
    }
    
    
    for(i in 1:perm){
      if(printRep & (i%%printRep==0)) cat("\n\nStarting optimization of the partiton",i,"of",perm,"partitions.\n")
      find.unique.par<-TRUE
      ununiqueParTested=0
      endFun<-FALSE
      while(find.unique.par){
        temppar<-part
        if(i>1){
          for(iPar in 1:length(temppar)){
            while(TRUE){
              tPar<-temppar[[iPar]]
              sel<-rbinom(n[iPar],size=1,prob = sharePerm)
              tPar[sel]<-sample(1:k[iPar],size=sum(sel),replace = TRUE)
              if(length(unique(tPar))==k[iPar]) break          
            }
            tPar->temppar[[iPar]]
          }
        } else {
          if(length(temppar)==1)temppar<-temppar[[1]]
          break
        }
        if(length(temppar)==1)temppar<-temppar[[1]]
        
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
  } # The else (muticore) option is not yet implemented
  # else {
  #   oneRep<-function(i,M,n,k,mingr,maxgr,addParam,rep, parGenFun,...){
  #     temppar<-parGenFun(n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam)
  #     #skip.par<-c(skip.par,list(temppar))
  #     
  #     tres <- try(stochBlock(M=M, clu=temppar,  ...))
  #     if(inherits(x = tres,what = "try-error")){
  #       tres<-list("try-error"=tres, err=Inf, startPart=temppar)
  #     }
  #     if(deleteMs){
  #       tres$M<-NULL
  #     }
  #     tres$best<-NULL
  #     return(list(tres))
  #   }
  #   
  #   if(!requireNamespace("doParallel")|!requireNamespace("doRNG")) useParLapply<-TRUE
  #   
  #   if(nCores==0){
  #     nCores<-detectCores()-1                    
  #   }
  #   
  #   pkgName<-utils::packageName()
  #   if(is.null(pkgName)) pkgName<-utils::packageName(environment(fun.by.blocks))
  #   if(useParLapply) {
  #     if(is.null(cl)) cl<-makeCluster(nCores)
  #     clusterSetRNGStream(cl)
  #     nC<-nCores
  #     #clusterExport(cl, varlist = c("kmBlock","kmBlockORP"))
  #     #clusterExport(cl, varlist = "kmBlock")
  #     clusterExport(cl, varlist = "pkgName", envir=environment())
  #     clusterEvalQ(cl, expr={requireNamespace(pkgName,character.only = TRUE)})
  #     res<-parLapplyLB(cl = cl,1:rep, fun = oneRep, M=M,n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam,rep=rep, parGenFun=parGenFun,...)
  #     if(stopcl) stopCluster(cl)
  #     res<-lapply(res,function(x)x[[1]])
  #   } else {
  #     requireNamespace("doParallel")
  #     requireNamespace("doRNG")
  #     if(!getDoParRegistered()|(getDoParWorkers()!=nCores)){
  #       if(!is.null(cl)) {
  #         #cl<-makeCluster(nCores)
  #         registerDoParallel(cl)
  #       } else registerDoParallel(nCores)
  #     }
  #     nC<-getDoParWorkers()
  #     
  #     res<-foreach(i=1:rep,.combine=c, .packages=pkgName) %dorng% oneRep(i=i,M=M,n=n,k=k,mingr=mingr,maxgr=maxgr,addParam=addParam,rep=rep, parGenFun=parGenFun,...)
  #     if(!is.null(cl) & stopcl) {
  #       registerDoSEQ()
  #       stopCluster(cl)
  #     }
  #   }
  #   err<-sapply(res,function(x)x$err)    
  # }
}

