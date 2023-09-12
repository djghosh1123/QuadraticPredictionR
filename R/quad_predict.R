
#' Calculate Linear or Quadratic Predictor
#' @param x Vector from GARCH process
#'
#' @param k.pred Number of future points to be predicted
#' @param type Type of Prediction, has to be quad or lin
#' @param past
#' @param iter Number of iterations
#'
#' @importFrom matrixcalc vech
#' @export

quad_predict = function(x,k.pred=1, past=20, type = "lin"){

  # improv = 0
  mu = mean(x)
  sd = sd(x)
  x = (x-mu)/sd

  sgm = NULL
  simx = x[(length(x)-past+1):(length(x))]
  sigma.xx = matrix(0,length(simx),length(simx))


  for(i in 1:(length(x)+k.pred-1)){
    v = general_cumulant(x,i-1)
    sgm = c(sgm,v)
  }

  sgm1 = sgm[1:past]

  for(i in 1:length(simx)){
    for(j in 1:length(simx)){
      if(j>=i){
        sigma.xx[i,j] = sgm1[abs(i-j)+1]
      }
      else{
        sigma.xx[i,j] = sigma.xx[j,i]
      }
    }
  }

  sigma.xy = matrix(0,length(simx),k.pred)
  for(i in 1:length(simx)){
    for(j in 1:k.pred){
      sigma.xy[i,j] = sgm[length(simx)-i+j+1]
    }
  }

  if(type=="lin"){
    y = t(sigma.xy)%*%solve(sigma.xx)%*%simx
  }

  if(type=="quad"){
    y = t(sigma.xy)%*%solve(sigma.xx)%*%simx

    ln = length(simx)
    cubic.mat = matrix(0,(ln*(ln+1))/2,ln)

    matrix.pair = NULL
    for(i in 1:length(x)){
      z1 = rep(i,length(x)-i+1)
      z2 = i:length(x)
      matrix.pair = rbind(matrix.pair, cbind(z1,z2))
    }

    for(i in 1:(ln*(ln+1)/2)){
      for(j in 1:ln){
        cubic.mat[i,j] = general_cumulant(x,as.numeric(c(
          matrix.pair[i,2] - matrix.pair[i,1],j- matrix.pair[i,1])))
      }
    }

    cubic.mat.y = matrix(0,ln*(ln+1)/2,k.pred)
    for(i in 1:(ln*(ln+1)/2)){
      for(j in 1:k.pred){
        cubic.mat.y[i,j] = general_cumulant(x,as.numeric(c(
          matrix.pair[i,2] - matrix.pair[i,1],length(simx)+j-matrix.pair[i,1])))
      }
    }

    quad.mat = matrix(0, ln*(ln+1)/2, ln*(ln+1)/2)

    for(i in 1:(ln*(ln+1)/2)){
      for(j in 1:(ln*(ln+1)/2)){
        if(j>=i){
          quad.mat[i,j] = general_cumulant(x,as.numeric(c(
            matrix.pair[i,2]- matrix.pair[i,1],
            matrix.pair[j,1]- matrix.pair[i,1],
            matrix.pair[j,2]- matrix.pair[i,1])))
        }
        else{
          quad.mat[i,j] = quad.mat[j,i]
        }
      }
    }

    schur.mat = quad.mat - cubic.mat%*%solve(sigma.xx)%*%
      t(cubic.mat)

    schur.mat = round(schur.mat,3)

    # e = eigen(schur.mat)
    # w = which(e$values<0)
    # e$values[w]=0.01
    # scmt = e$vectors%*%diag(e$values)%*%t(e$vectors)



    M = rbind(cbind(sigma.xx,t(cubic.mat)),cbind(cubic.mat,quad.mat))

    sigma.vech.e = cubic.mat.y - cubic.mat%*%solve(sigma.xx)%*%
      (sigma.xy)

    cmv = cubic.mat%*%(solve(sigma.xx)%*%sigma.xy)

    bet1 = solve(sigma.xx)%*%sigma.xy - (solve(sigma.xx)%*%t(cubic.mat))%*%
      solve(schur.mat)%*%(cubic.mat.y - cmv)

    bet2 = solve(schur.mat)%*%(cubic.mat.y - cmv)

    bet = rbind(bet1,bet2)

    x.new = rbind(t(t(simx)), t(t(matrixcalc::vech(simx%*%t(simx)))))

    #y = y - t(sigma.vech.e)%*%solve(schur.mat)%*%
    #   (vech(x%*%t(x))- cubic.mat%*%solve(sigma.xx)%*%x)

    y = t(bet)%*%x.new/10

    # improv = t(sigma.vech.e)%*%solve(schur.mat)%*%sigma.vech.e

  }

  return(y*sd + mu)

  # return(list("pred"=y,"improv" = improv))
}

