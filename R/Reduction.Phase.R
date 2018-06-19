
######################################################
### Reduction.Phase : Return Matrix.Selection and List.Variable.Selection
######################################################

Reduction.Phase = function(X,Y,family=gaussian,log.transf=FALSE, dmHC=NULL,vector.signif=NULL,seed.HC = NULL){

  n = nrow(X)
  d = ncol(X)

  ### Hypercube dimension Not specified by the user
  if(is.null(dmHC)){
    v=2:5
    dmHC = v[ceiling(d^(1/v)) <= 15 & ceiling(d^(1/v))>=10]
    if(dmHC==0){
      stop('Hypercube dimension was not specified and conditioned to the number of variables in the design matrix could not the calculated.','\n',
           'Consider delete/include/transform some variables...')
    }
  }

  ### Error messages
  if(dmHC>5){
    stop('Sorry, this version only support cube with maximum 5 dimensions! More dimensions will be avaliable in the next version...')
  }

  if(!is.null(dmHC) & !is.null(vector.signif) & (length(vector.signif)+1)!=dmHC){
    stop(paste('The lenght of vector.signif is not the same as the HC dimension',dmHC,'.'))
  }

  ### Design Matrix Transformation
  if(log.transf==TRUE){
    lX = log(X)
    mLX = colMeans(lX)
    X = scale(lX, center=TRUE, scale=FALSE)
  }

  ### Signif vectors
  if(is.null(vector.signif)){
    signif.Default =TRUE
  } else{
    signif.Default =FALSE
    vector.signif = c(NA,vector.signif)
  }

  highest.dmHC = dmHC

  ### Outputs
  Matrix.Selection = list()
  List.Selection = list()

  aux.dmHC5 = aux.dmHC4 = aux.dmHC3 = aux.dmHC2 = 'N'



  ########## case in which dmHC=5 ##########
  if(dmHC==5){

    if(signif.Default==TRUE & highest.dmHC==5){
      aux.signif = 2
    } else{
      aux.signif = vector.signif[5]
    }

    dimHC5 = ceiling(d^(1/5))
    nearestHypercube5=dimHC5^5
    remainderHC5=nearestHypercube5-d

    if(!is.null(seed.HC)){
      set.seed(seed.HC)
      hypercube5=array(sample(c((1:d),rep(0,remainderHC5))),dim=c(dimHC5,dimHC5,dimHC5,dimHC5,dimHC5))
    } else{
      hypercube5=array(c((1:d),rep(0,remainderHC5)),dim=c(dimHC5,dimHC5,dimHC5,dimHC5,dimHC5))
    }
    ### intermediate step to handle 0`s here
    hypercubeSelect=array(0,dim=c(dimHC5,dimHC5,dimHC5,dimHC5,dimHC5))

    for(ind5 in 1:dim(hypercube5)[5]){
      for(ind4 in 1:dim(hypercube5)[4]){
        for(ind3 in 1:dim(hypercube5)[3]){
          for(indR in 1:dim(hypercube5)[1]){
            if(length(which(hypercube5[indR,,ind3,ind4,ind5]!=0))>0){
              idx.aux = which(hypercube5[indR,,ind3,ind4,ind5]!=0)
              subsetX = cbind(X[,hypercube5[indR,idx.aux,ind3,ind4,ind5]])
              idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif)

              if(length(idxSelected$idxSelected)>0){
                hypercubeSelect[indR,idxSelected$idxSelected,ind3,ind4,ind5]=hypercubeSelect[indR,idxSelected$idxSelected,ind3,ind4,ind5]+1
              }
            }
          } # indR
          for(indC in 1:dim(hypercube5)[2]){
            if(length(which(hypercube5[,indC,ind3,ind4,ind5]!=0))>0){
              idx.aux = which(hypercube5[,indC,ind3,ind4,ind5]!=0)
              subsetX = cbind(X[,hypercube5[idx.aux,indC,ind3,ind4,ind5]])
              idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif)

              if(length(idxSelected$idxSelected)>0){
                hypercubeSelect[idxSelected$idxSelected,indC,ind3,ind4,ind5]=hypercubeSelect[idxSelected$idxSelected,indC,ind3,ind4,ind5]+1
              }
            }
          } # indC
        } # ind3
        for(indR in 1:dim(hypercube5)[1]){
          for(indC in 1:dim(hypercube5)[2]){
            if(length(which(hypercube5[indR,indC,,ind4,ind5]!=0))>0){
              idx.aux = which(hypercube5[indR,indC,,ind4,ind5]!=0)
              subsetX = cbind(X[,hypercube5[indR,indC,idx.aux,ind4,ind5]])
              idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif)

              if(length(idxSelected$idxSelected)>0){
                hypercubeSelect[indR,indC,idxSelected$idxSelected,ind4,ind5]=hypercubeSelect[indR,indC,idxSelected$idxSelected,ind4,ind5]+1
              }
            }
          } # indR
        } # indC
      } # ind4
      for(ind3 in 1:dim(hypercube5)[3]){
        for(indR in 1:dim(hypercube5)[1]){
          for(indC in 1:dim(hypercube5)[2]){
            if(length(which(hypercube5[indR,indC,ind3,,ind5]!=0))>0){
              idx.aux = which(hypercube5[indR,indC,ind3,,ind5]!=0)
              subsetX = cbind(X[,hypercube5[indR,indC,ind3,idx.aux,ind5]])
              idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif)

              if(length(idxSelected$idxSelected)>0){
                hypercubeSelect[indR,indC,ind3,idxSelected$idxSelected,ind5]=hypercubeSelect[indR,indC,ind3,idxSelected$idxSelected,ind5]+1
              }
            }
          } # indC
        } # indR
      } # ind3
    } # ind5
    # now traverse in the 5th dimension
    for(ind4 in 1:dim(hypercube5)[4]){
      for(ind3 in 1:dim(hypercube5)[3]){
        for(indR in 1:dim(hypercube5)[1]){
          for(indC in 1:dim(hypercube5)[2]){
            if(length(which(hypercube5[indR,indC,ind3,ind4,]!=0))>0){
              idx.aux = which(hypercube5[indR,indC,ind3,ind4,]!=0)
              subsetX = cbind(X[,hypercube5[indR,indC,ind3,ind4,idx.aux]])
              idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif)

              if(length(idxSelected$idxSelected)>0){
                hypercube5[indR,indC,ind3,ind4,idxSelected$idxSelected]=hypercube5[indR,indC,ind3,ind4,idxSelected$idxSelected]+1
              }
            }
          } # indC
        } # indR
      } # ind3
    }

    setSelected4Times5 = hypercube5[which(hypercubeSelect>3, arr.ind = TRUE)]
    setSelected3Times5 = hypercube5[which(hypercubeSelect>2, arr.ind = TRUE)]
    setSelected2Times5 = hypercube5[which(hypercubeSelect>1, arr.ind = TRUE)]

    numSelected4times5=length(setSelected4Times5)
    numSelected3times5=length(setSelected3Times5)
    numSelected2times5=length(setSelected2Times5)

    Matrix.Selection[[paste('Hypercube with dim',dmHC)]] = c(numSelected2times5,numSelected3times5,numSelected4times5)
    names(Matrix.Selection[[paste('Hypercube with dim',dmHC)]]) = c('numSelected2','numSelected3','numSelected4')

    List.Selection[[paste('Hypercube with dim',dmHC)]][[paste('numSelected2')]] = setSelected2Times5
    List.Selection[[paste('Hypercube with dim',dmHC)]][[paste('numSelected3')]] = setSelected3Times5
    List.Selection[[paste('Hypercube with dim',dmHC)]][[paste('numSelected4')]] = setSelected4Times5

    aux.dmHC5 <- 'Y' #readline(cat("Reduction of dimension 5 done!", "\n", length(setSelected4Times5),
    #      "Variables selected at least 4 times","\n",length(setSelected3Times5),
    #     "Variables selected at least 3 times","\n",length(setSelected2Times5),
    #    "Variables selected at least 2 times","\n", "Wanna proceed with reduction?[Y/N]"))

  }

  ########## case in which dmHC=4 ##########
  if(dmHC==4 | aux.dmHC5=='Y'){

    if(signif.Default==TRUE & highest.dmHC==4){
      aux.signif = 2
    } else if(signif.Default==TRUE & highest.dmHC>4){
      aux.signif = 0.01
    } else{
      aux.signif = vector.signif[4]
    }

    dimHC = ceiling(d^(1/4))
    nearestHypercube=dimHC^4
    remainderHC=nearestHypercube-d
    if(!is.null(seed.HC)){
      set.seed(seed.HC)
      hypercube=array(sample(c((1:d),rep(0,remainderHC))),dim=c(dimHC,dimHC,dimHC,dimHC))
    } else{
      hypercube=array(c((1:d),rep(0,remainderHC)),dim=c(dimHC,dimHC,dimHC,dimHC))
    }
    ### intermediate step to handle 0`s here
    hypercubeSelect=array(0,dim=c(dimHC,dimHC,dimHC,dimHC))

    for(ind4 in 1:dim(hypercube)[4]){
      for(ind3 in 1:dim(hypercube)[3]){
        for(indR in 1:dim(hypercube)[1]){
          if(length(which(hypercube[indR,,ind3,ind4]!=0))>0){
            idx.aux = which(hypercube[indR,,ind3,ind4]!=0)
            subsetX = cbind(X[,hypercube[indR,idx.aux,ind3,ind4]])
            idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif)

            if(length(idxSelected$idxSelected)>0){
              hypercubeSelect[indR,idxSelected$idxSelected,ind3,ind4]=hypercubeSelect[indR,idxSelected$idxSelected,ind3,ind4]+1
            }
          }
        } # indR
        for(indC in 1:dim(hypercube)[2]){
          if(length(which(hypercube[,indC,ind3,ind4]!=0))>0){
            idx.aux = which(hypercube[,indC,ind3,ind4]!=0)
            subsetX = cbind(X[,hypercube[idx.aux,indC,ind3,ind4]])
            idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif)

            if(length(idxSelected$idxSelected)>0){
              hypercubeSelect[idxSelected$idxSelected,indC,ind3,ind4]=hypercubeSelect[idxSelected$idxSelected,indC,ind3,ind4]+1
            }
          }
        } # indC
      } # ind3
      for(indR in 1:dim(hypercube)[1]){
        for(indC in 1:dim(hypercube)[2]){
          if(length(which(hypercube[indR,indC,,ind4]!=0))>0){
            idx.aux = which(hypercube[indR,indC,,ind4]!=0)
            subsetX = cbind(X[,hypercube[indR,indC,idx.aux,ind4]])
            idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif)

            if(length(idxSelected$idxSelected)>0){
              hypercubeSelect[indR,indC,idxSelected$idxSelected,ind4]=hypercubeSelect[indR,indC,idxSelected$idxSelected,ind4]+1
            }
          }
        } # indC
      } # indR
    } # ind4
    # now traverse in the 4th dimension
    for(ind3 in 1:dim(hypercube)[3]){
      for(indR in 1:dim(hypercube)[1]){
        for(indC in 1:dim(hypercube)[2]){
          if(length(which(hypercube[indR,indC,ind3,]!=0))>0){
            idx.aux = which(hypercube[indR,indC,ind3,]!=0)
            subsetX = cbind(X[,hypercube[indR,indC,ind3,idx.aux]])
            idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif)

            if(length(idxSelected$idxSelected)>0){
              hypercubeSelect[indR,indC,ind3,idxSelected$idxSelected]=hypercubeSelect[indR,indC,ind3,idxSelected$idxSelected]+1
            }
          }
        } # indC
      } # indR
    } # ind3

    setSelected4Times = hypercube[which(hypercubeSelect>3, arr.ind = TRUE)]
    setSelected3Times = hypercube[which(hypercubeSelect>2, arr.ind = TRUE)]
    setSelected2Times = hypercube[which(hypercubeSelect>1, arr.ind = TRUE)]

    numSelected4times5=length(setSelected4Times)
    numSelected3times5=length(setSelected3Times)
    numSelected2times5=length(setSelected2Times)

    Matrix.Selection[[paste('Hypercube with dim',4)]] = c(numSelected2times5,numSelected3times5,numSelected4times5)
    names(Matrix.Selection[[paste('Hypercube with dim',4)]]) = c('numSelected2','numSelected3','numSelected4')

    List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected2')]] = setSelected2Times
    List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected3')]] = setSelected3Times
    List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected4')]] = setSelected4Times


    aux.dmHC4 <- 'Y' #readline(cat("Reduction of dimension 4 done!", "\n", length(setSelected4Times5),
    #        "Variables selected at least 4 times","\n",length(setSelected3Times5),
    #        "Variables selected at least 3 times","\n",length(setSelected2Times5),
    #       "Variables selected at least 2 times","\n", "Wanna proceed with reduction?[Y/N]"))
  }

  ########## case in which dmHC=3 ##########
  if(dmHC==3 | aux.dmHC4=='Y'){

    if(signif.Default==TRUE & highest.dmHC==3){
      aux.signif = 2
    } else if(signif.Default==TRUE & highest.dmHC>3){
      aux.signif = 0.01
    } else{
      aux.signif = vector.signif[3]
    }

    dimHC = ceiling(d^(1/3))
    nearestHypercube=dimHC^3
    remainderHC=nearestHypercube-d
    if(!is.null(seed.HC)){
      set.seed(seed.HC)
      hypercube=array(sample(c((1:d),rep(0,remainderHC))),dim=c(dimHC,dimHC,dimHC))
    } else{
      hypercube=array(c((1:d),rep(0,remainderHC)),dim=c(dimHC,dimHC,dimHC))
    }
    ### intermediate step to handle 0`s here
    hypercubeSelect=array(0,dim=c(dimHC,dimHC,dimHC))

    for(indL in 1:dim(hypercube)[3]){
      for(indR in 1:dim(hypercube)[1]){
        if(length(which(hypercube[indR,,indL]!=0))>0){
          idx.aux = which(hypercube[indR,,indL]!=0)
          subsetX = cbind(X[,hypercube[indR,idx.aux,indL]])
          idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif)

          if(length(idxSelected$idxSelected)>0){
            hypercubeSelect[indR,idxSelected$idxSelected,indL]=hypercubeSelect[indR,idxSelected$idxSelected,indL]+1
          }
        }
      } # indR
      for(indC in 1:dim(hypercube)[2]){
        if(length(which(hypercube[,indC,indL]!=0))>0){
          idx.aux = which(hypercube[,indC,indL]!=0)
          subsetX = cbind(X[,hypercube[idx.aux,indC,indL]])
          idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif)

          if(length(idxSelected$idxSelected)>0){
            hypercubeSelect[idxSelected$idxSelected,indC,indL]=hypercubeSelect[idxSelected$idxSelected,indC,indL]+1
          }
        }
      } # indC
    } # indL
    for(indR in 1:dim(hypercube)[1]){
      for(indC in 1:dim(hypercube)[2]){
        if(length(which(hypercube[indR,indC,]!=0))>0){
          idx.aux = which(hypercube[indR,indC,]!=0)
          subsetX = cbind(X[,hypercube[indR,indC,idx.aux]])
          idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif)

          if(length(idxSelected$idxSelected)>0){
            hypercubeSelect[indR,indC,idxSelected$idxSelected]=hypercubeSelect[indR,indC,idxSelected$idxSelected]+1
          }
        }
      } # indC
    } # indR

    setSelected2Times = hypercube[which(hypercubeSelect>1, arr.ind = TRUE)]
    setSelected1Times = hypercube[which(hypercubeSelect>0, arr.ind = TRUE)]

    numSelected2times=length(setSelected2Times)
    numSelected1times=length(setSelected1Times)

    Matrix.Selection[[paste('Hypercube with dim',3)]] = c(numSelected1times,numSelected2times)
    names(Matrix.Selection[[paste('Hypercube with dim',3)]]) = c('numSelected1','numSelected2')

    List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected1')]] = setSelected1Times
    List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected2')]] = setSelected2Times

    aux.dmHC3 <- 'Y' #readline(cat("Reduction of dimension 3 done!", "\n", length(setSelected2Times),
    #         "Variables selected at least 2 times","\n",length(setSelected1Times),
    #        "Variables selected at least 1 times","\n", "Wanna proceed with reduction?[Y/N]"))
  }

  ########## case in which dmHC=2 ##########
  if(dmHC==2 | aux.dmHC3=='Y'){

    if(signif.Default==TRUE & highest.dmHC==2){
      aux.signif = 2
    } else if(signif.Default==TRUE & highest.dmHC>2){
      aux.signif = 0.01
    } else{
      aux.signif = vector.signif[2]
    }

    dimHC = ceiling(d^(1/2))
    nearestHypercube=dimHC^2
    remainderHC=nearestHypercube-d
    if(!is.null(seed.HC)){
      set.seed(seed.HC)
      hypercube=array(sample(c((1:d),rep(0,remainderHC))),dim=c(dimHC,dimHC))
    } else{
      hypercube=array(c((1:d),rep(0,remainderHC)),dim=c(dimHC,dimHC))
    }
    ### intermediate step to handle 0`s here
    hypercubeSelect=array(0,dim=c(dimHC,dimHC))

    for(indR in 1:dim(hypercube)[1]){
      if(length(which(hypercube[indR,]!=0))>0){
        idx.aux = which(hypercube[indR,]!=0)
        subsetX = cbind(X[,hypercube[indR,idx.aux]])
        idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif)

        if(length(idxSelected$idxSelected)>0){
          hypercubeSelect[indR,idxSelected$idxSelected]=hypercubeSelect[indR,idxSelected$idxSelected]+1
        }
      }
    } # indR
    for(indC in 1:dim(hypercube)[2]){
      if(length(which(hypercube[,indC]!=0))>0){
        idx.aux = which(hypercube[,indC]!=0)
        subsetX = cbind(X[,hypercube[idx.aux,indC]])
        idxSelected = glmRoutine(subsetX,Y,family=family,significance = aux.signif)

        if(length(idxSelected$idxSelected)>0){
          hypercubeSelect[idxSelected$idxSelected,indC]=hypercubeSelect[idxSelected$idxSelected,indC]+1
        }
      }
    } # indC

    setSelected2Times = hypercube[which(hypercubeSelect>1, arr.ind = TRUE)]
    setSelected1Times = hypercube[which(hypercubeSelect>0, arr.ind = TRUE)]

    numSelected2times=length(setSelected2Times)
    numSelected1times=length(setSelected1Times)

    Matrix.Selection[[paste('Hypercube with dim',2)]] = c(numSelected1times,numSelected2times)
    names(Matrix.Selection[[paste('Hypercube with dim',2)]]) = c('numSelected1','numSelected2')

    List.Selection[[paste('Hypercube with dim',2)]][[paste('numSelected1')]] = setSelected1Times
    List.Selection[[paste('Hypercube with dim',2)]][[paste('numSelected2')]] = setSelected2Times

    aux.dmHC2 <- 'Y' #readline(cat("Reduction of dimension 2 done!", "\n", length(setSelected2Times),
    #"Variables selected at least 2 times","\n",length(setSelected1Times),
    #"Variables selected at least 1 times","\n", "Wanna proceed with reduction?[Y/N]"))
  }

  return(list("Matrix.Selection" = Matrix.Selection, "List.Selection" = List.Selection))
}
