######################################################
### Exploratory.Phase : Return the significant Squared terms and iteration terms
######################################################

Exploratory.Phase = function(X, Y, list.reduction, family=gaussian, signif=0.05, silent=TRUE){

  n = nrow(X)

  #list.reduction  = aux.XX2$List.Selection$`Hypercube with dim 3`$numSelected2
  idx.combinations = t(combn(list.reduction,2))
  mat.select.SQ = mat.select.INTER = NULL

  for(ii in 1:length(list.reduction)){
    ### squared
    xTmp = matrix(0,n,2)
    xTmp[,1] = X[,list.reduction[ii]]
    xTmp[,2] = X[,list.reduction[ii]]^2
    stats = glm(Y~xTmp, family = family)
    pVals = coef(summary(stats))[,4]
    pVals = pVals[-1]
    tSt = coef(summary(stats))[,3]
    tSt = tSt[-1]
    if(pVals[2]<signif){
      mat.select.SQ = c(mat.select.SQ,list.reduction[ii])
    }
  }

  for(ii in 1:nrow(idx.combinations)){
    ### iteraction
    xTmp = matrix(0,n,3)
    xTmp[,1] = X[,idx.combinations[ii,1]]
    xTmp[,2] = X[,idx.combinations[ii,2]]
    xTmp[,3] = X[,idx.combinations[ii,1]] * X[,idx.combinations[ii,2]]
    stats = glm(Y~xTmp, family = family)
    pVals = coef(summary(stats))[,4]
    pVals = pVals[-1]
    tSt = coef(summary(stats))[,3]
    tSt = tSt[-1]
    if(pVals[3]<signif){
      mat.select.INTER = rbind(mat.select.INTER,c(idx.combinations[ii,1],idx.combinations[ii,2]))
    }
  }

  if(silent==FALSE){

    if(is.null(mat.select.INTER)){
      stop('No variables selected with interaction! Please increase the significance level.')
    }

    ## deciding upon type of Y
    type.var = ifelse(length(unique(Y))==2,"B","C")

    mat.response.INTER = NULL
    for(i in 1:nrow(mat.select.INTER)){

      data.res = data.frame("X1"=X[,mat.select.INTER[i,1]],
                            "X2"=X[,mat.select.INTER[i,2]], "Y"= Y)
      names(data.res)[1] = mat.select.INTER[i,1]
      names(data.res)[2] = mat.select.INTER[i,2]

      if(type.var =="C"){
        ## Continuous response
        print(ggplot(data.res, aes(x=data.res[,1], y=data.res[,2], color=Y)) + geom_point() +
                labs(x = mat.select.INTER[i,1], y = mat.select.INTER[i,2]) + scale_color_gradientn(colours = rainbow(2)))
      }

      if(type.var =="B"){
        ## Continuous response
        respost = unique(Y)
        print(ggplot(data.res, aes(x=data.res[,1], y=data.res[,2], color=factor(Y))) + geom_point() +
          lims(colour = c(toString(respost[1]), toString(respost[2]))) + labs(x = mat.select.INTER[i,1], y = mat.select.INTER[i,2], colour = "Y"))

      }

      controle=0
      while (controle!=1) {
        answ = readline(prompt = "Would you like to discard this interaction term? [Y/N] > ")
        if(answ=="Y" | answ=="N"){
          controle=1
        }
      }
      mat.response.INTER[i] = answ
    }

    mat.select.INTER = mat.select.INTER[which(mat.response.INTER=="N"),]

  }


  return(list("mat.select.SQ"=mat.select.SQ,"mat.select.INTER"=mat.select.INTER))

}
