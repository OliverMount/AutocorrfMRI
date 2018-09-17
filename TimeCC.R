TimeCC<-function(x,y)
{
  # Time Cross correlation function  
  # Asymtotically unbiased estimator
  # x and y are assumed to be of same length
  # Cross correlation is not symmetric (careful to deal with the negative lags)
  
  # Initializations
  M=length(x)
  rp <- array(0, dim = c(M))   # Positive part
  rn <- array(0, dim = c(M-1)) # Negative part
  
  # For positive lags
  for (i in 1:M)
  {
    rp[i]=(1/(M))*(x[i:M]%*%y[1:(M-i+1)])   # dot product in R
  }
  
  # For negative lags
  for (i in 1:M)
  {
    rn[i]=(1/(M))*(y[i:M]%*%x[1:(M-i+1)])   # dot product in R
  }
  
  # Combined negative and positive parts returned
  return(c(rev(rn[2:M]),rp))
}