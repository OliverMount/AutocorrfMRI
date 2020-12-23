VarSPCC<- function(R1,R2)
{

  # R1 and R2 are two time series
  N= length(R1)
  # Initializations
  D<-array(0,dim=c(3,3))
 
  # Estimate the time series correlation
  rx=TimeAC(R1)
  ry=TimeAC(R2)

 #Estimate the time series cross-correlation
  rxy=TimeCC(R1,R2) 
  ryx=TimeCC(R2,R1)
 

  rx0=rx[N]  # rx(0)
  ry0=ry[N]  # ry(0)
  rxy0=rxy[N]  # rxy(0)
 
 
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
 
  # Variance of the SPCC
  va=(colSums(C* (D %*% C)))/N
 
  return(va)   # Variance of the SPCC   
}
