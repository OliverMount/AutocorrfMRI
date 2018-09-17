DOFofAR1 <- function(Ch,alpha,beta,N,rho,varv,varw)
{
  
  # This function returns the variance of the SPC for two AR(1) 
  # correlated time series with parameters alpha and beta
  
  # Initializations
  Cmtx<-array(0,dim=c(3,3))
  
  # Input side
  rv=var1*(alpha^(abs(seq(-(N-1),(N-1),by=1)))/(1-(alpha^2))) # Impulse ACF
  rw=var2*(beta^(abs(seq(-(N-1),(N-1),by=1)))/(1-(beta^2)))
  
  # rvw
  n=0:(N-1) 
  right=((1-((alpha*beta)^(N-n)))/(1-(alpha*beta)))*(alpha^n)
  left=((1-((alpha*beta)^(N-(rev(n)))))/(1-(alpha*beta)))*(beta^(rev(n)))
  rvw=rho*sqrt(var1)*sqrt(var2)*c(left[1:N-1],right)
  
  # rwv
  right=((1-((alpha*beta)^(N-n)))/(1-(alpha*beta)))*(beta^n)
  left=((1-((alpha*beta)^(N-(rev(n)))))/(1-(alpha*beta)))*(alpha^(rev(n)))
  rwv=rho*sqrt(var1)*sqrt(var2)*c(left[1:N-1],right)
  
  # Output side
  rx=conv(rv,Ch)
  ry=conv(rw,Ch)
  rxy=conv(Ch,rvw)  # output cross correlation calculation
  ryx=conv(Ch,rwv)  # output cross correlation calculation
  
  
  C=c(-rvw[N]/(2*(sqrt((rv[N]^3)*rw[N]))),
      -rvw[N]/(2*(sqrt(rv[N]*(rw[N]^3)))),
      1/sqrt(rv[N]*rw[N]))
  
  D=array(0,dim=c(3,3))
  # Formation of Covariance matrix
  # Factor 2 is removed on the diagonals to have the correct transpose
  D[1,1]=sum(rv^2)
  D[1,2]=2*sum(rvw^2)
  D[1,3]=2*sum(rwv*rv)
  D[2,2]=sum(rw^2)
  D[2,3]=2*sum(rvw*rw)
  D[3,3]=0.5*(sum(rv*rw)+sum(rvw*rwv))
  D=(D+t(D)) 
  
  va=(colSums(C* (D %*% C)))/N
  me=(rho*sqrt(1-(alpha^2))*sqrt(1-(beta^2)))/(1-(alpha*beta))
  
  return(c(me,va))
}



