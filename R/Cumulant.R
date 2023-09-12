

#' Calculate Second order Autocumulant
#'
#' @param lag The lags at which to calculate cumulants, a scalar
#' @param x Sample from which autocumulant to be calculated
#' @export


second_order_cumulant = function(x,lag){

  y = x[(abs(lag)+1):length(x)]
  x.new = x[1:length(y)]
  mnt = mean(y*x.new)

  return(mnt)
}


#' Calculate Third order Autocumulant
#'
#' @param x Sample from which autocumulant to be calculated
#' @param lag The lags at which to calculate cumulants, a vector of length 2
#' @export


third_order_cumulant = function(x,lag){
  #x = garch.sim(alpha,n,b,g,...)
  if(lag[1]<0 || lag[2]<0){
    z1 = min(lag)
    z2 = max(lag)
    lag = c(-z1, z2-z1)
  }

  y1 = x[(min(lag)+1):length(x)]
  y2 = x[(max(lag)+1):length(x)]
  x.new = x[1:length(y2)]
  y1 = y1[1:length(y2)]
  mnt = mean(y1*x.new*y2)


  return(mnt)
}


#' Calculate Fourth order Autocumulant
#'
#' @param x Sample
#' @param lag The lags at which to calculate cumulants, a vector of length 3
#'
#' @export


fourth_order_cumulant = function(x,lag){
  #x = garch.sim(alpha,n,b,g,...)
  if(lag[1]<0 || lag[2]<0 || lag[3]<0){
    z = sort(lag)
    lag = c(-z[1], z[2]-z[1], z[3]-z[1])
  }

  lag = sort(lag)
  y1 = x[(lag[1]+1):length(x)]
  y2 = x[(lag[2]+1):length(x)]
  y3 = x[(lag[3]+1):length(x)]
  x.new = x[1:length(y3)]
  y1 = y1[1:length(y3)]
  y2 = y2[1:length(y3)]
  mnt = mean(y1*x.new*y2*y3)

  mnt = mnt + second_order_cumulant(x,lag[1])*
    second_order_cumulant(x,lag[2]-lag[3]) +
    second_order_cumulant(x,lag[2])*
    second_order_cumulant(x,lag[3]-lag[1]) +
    second_order_cumulant(x,lag[3])*
    second_order_cumulant(x,lag[1]-lag[2])


  return(mnt)
}

#' Autocumulant Calculation (second, third or fourth)
#' @param x Sample
#' @param lag Lag vector
#' @param iter Number of Monte Carlo Simulations
#' @export


general_cumulant = function(x,lag){

    if(length(lag)==1)
      mnt = second_order_cumulant(x,lag)
    else if(length(lag)==2)
      mnt = third_order_cumulant(x,lag)
    else if(length(lag)==3)
      mnt = fourth_order_cumulant(x,lag)
    else
      stop("Maximum length of lag can be 3")

    return(mnt)
}




