######################################################
### glmRoutine
######################################################

glmRoutine = function(subsetX,Y,family=gaussian,intercept=TRUE,significance=NULL){

  n = nrow(subsetX)

  if(intercept==TRUE){
    intercept=1
  } else{
    intercept=0
  }
  count = 0

  if(intercept==0){
    stats = glm(Y~subsetX-1, family = family)
    pVals = coef(summary(stats))[,4]
  } else{
    stats = glm(Y~subsetX, family = family)
    pVals = coef(summary(stats))[,4]
    pVals = pVals[-1]
  }

  if(length(which(pVals>0.99))==length(pVals) | length(which(pVals==0))==length(pVals)){
    sTmp=ncol(subsetX)
    sB1=floor(sTmp/2)
    block1=subsetX[,1:sB1]
    block2=subsetX[,(sB1+1):ncol(subsetX)]
    stats = glm(Y~block1, family = family)
    pVals1 = coef(summary(stats))[,4]
    pVals1 = pVals1[-1]
    if(length(which(pVals1>0.99))==length(pVals1)){
      pVals1 = 0.0001*rep(1,length(pVals1))
      count = count+1
    }
    stats = glm(Y~block2, family = family)
    pVals2=coef(summary(stats))[,4]
    pVals2 = pVals2[-1]
    if(length(which(pVals2>0.99))==length(pVals2)){
      pVals2 = 0.0001*rep(1,length(pVals2))
      count = count+1
    }
    pVals=c(pVals1,pVals2)
  }

  if(significance==2){
    idxSelected=which(pVals %in% sort(pVals)[1:2])
  } else if(significance==1){
    idxSelected=which(pVals %in% sort(pVals)[1])
  }  else{
    idxSelected=which(pVals<=significance)
  }

  return(list("idxSelected"=idxSelected,"count"=count))
}

