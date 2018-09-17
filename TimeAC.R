TimeAC <- function(x)
{

  M <-length(x)
  r <- array(0, dim = c(M))  # One-sided autocorrelation

  for (i in 1:M) 
  {
    r[i]=(x[1:(M-i+1)]%*%x[i:M])    # Dot product in r
  }
  
  rr=c(rev(r[2:M]),r[1:M]) 
  return(rr/M)

}