######################################################
### Model.Selection.Phase : Return the model confidence set
######################################################

Model.Selection.Phase = function(X,Y, list.reduction, family=binomial, signif=0.05, square.terms=NULL, interaction.terms=NULL, modelSize=NULL){

  if(!is.null(square.terms)){
    X.SQ = X[,square.terms]^2
    SQ.names = paste(square.terms,"^2")
  }

  if(is.null(square.terms)){
    X.SQ = SQ.names = NULL
  }

  if(!is.null(interaction.terms)){
    X.ITER = iter.names = NULL
    for(ii in 1:nrow(interaction.terms)){
      X.ITER = cbind(X.ITER,X[,interaction.terms[ii,1]]*X[,interaction.terms[ii,2]])
      iter.names = c(iter.names,paste(interaction.terms[ii,1],"-",interaction.terms[ii,2]))
    }
  }

  if(is.null(interaction.terms)){
    X.ITER = iter.names = NULL
  }

  #### como fica daqui para baixo com os termos NULL???

  X.full = cbind(X[,list.reduction],X.SQ,X.ITER)
  setSelected = c(list.reduction,SQ.names,iter.names)

  s=ncol(X.full)
  if(is.null(modelSize)){ modelSize=min(5,length(setSelected)) } ## check this with Heather

  ### Error message
  if(modelSize>7){
    stop('Sorry, this version only support model sizes<8')
  }

  n=nrow(X.full)
  residBig = residuals(glm(Y~X.full, family = family))
  RSSBig = t(residBig) %*%residBig
  dfBig = n-(length(setSelected)+1)

  goodModels = list()

  for (j in 1:modelSize){
    combinationMatrix=t(combn(1:length(setSelected),j))
    combinationMatrixNames=cbind(t(combn(setSelected,j)))
    logicFitVectorF=matrix(0,nrow(combinationMatrix))
    for(l in 1:nrow(combinationMatrix)){
      XSelect = X.full[,combinationMatrix[l,]]
      residSmall = residuals(glm(Y~XSelect, family = family))
      RSSSmall = t(residSmall) %*% residSmall
      diffRSS = RSSSmall - RSSBig
      dfSmall = n-(ncol(cbind(XSelect))+1)
      fCrit = qf(1-signif,df1 = (dfSmall-dfBig),df2 = dfBig)
      if(is.nan(fCrit)){
        logicFitVectorF[l]=1
      } else if((diffRSS/(dfSmall-dfBig))/(RSSBig/dfBig)<=fCrit){
        logicFitVectorF[l]=1
      }
    }
    if(j==1){goodModels[[paste('Model Size', j)]] = as.matrix(combinationMatrixNames[which(logicFitVectorF>0),])
    } else{
      goodModels[[paste('Model Size', j)]] = combinationMatrixNames[which(logicFitVectorF>0),]
    }
  }

  return(list("goodModels"=goodModels))

}
