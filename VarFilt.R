VarFilt <- function(Ch,alpha,beta,N,rho,varv,varw)
{
  
  # Returns the variance of SPC of AR(1) model including the filter 
  
  # Initializations
  CV=(2*N)-2    # Center value
  D<-array(0,dim=c(3,3))
  
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
  #rx0=rx[CV]
  rx0=rx[ceiling(length(rx)/2)]
  ry=conv(rw,Ch)
  #ry0=ry[CV]
  ry0=ry[ceiling(length(ry)/2)]
  
  rxy=conv(Ch,rvw)  # output cross correlation calculation
  #rxy0=rxy[CV]
  rxy0=rxy[ceiling(length(rxy)/2)]
  ryx=conv(Ch,rwv)  # output cross correlation calculation
  
  
  # Formation of Covariance matrix
  
  C=c(-rxy0/(2*(sqrt((rx0^3)*ry0))),
      -rxy0/(2*(sqrt(rx0*(ry0^3)))),
      1/sqrt(rx0*ry0))
  
  
  # Formation of Covariance matrix
  # Factor 2 is removed on the diagonals to have the correct transpose
  D[1,1]=sum(rx^2)
  D[1,2]=2*sum(rxy^2)
  D[1,3]=2*sum(rxy*rx)
  D[2,2]=sum(ry^2)
  D[2,3]=2*sum(rxy*ry)
  D[3,3]=0.5*(sum(rx*ry)+sum(rxy*ryx))
  D=(D+t(D)) 
  
  va=(colSums(C* (D %*% C)))/N
  me=rxy0/sqrt(rx0*ry0)
  
  return(c(me,va))   
}
