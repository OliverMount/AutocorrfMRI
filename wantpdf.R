wantpdf <- function(sig2esti,Ncells)
{
  # Empirical PDF estimator
  #Function WANTPDF Calculates the PDF of the signal sig2esti
  #N cells Specifies the number of cells for Histogram estimation
  # and if theoretical PDF (Gaussian) is needed, input the mean and variance 
  #function [p, histo, ptheory]=wantpdf(sig2esti,Ncells,m,variance)
  
  # Finding the PDF
  M=length(sig2esti)
  a=max(sig2esti)
  b=min(sig2esti)
  range=a-b
  Hepsi=(range/(2*Ncells))
  values=seq((b+Hepsi),a,by=(2*Hepsi))
  
  histo<-p <-array(0,dim=c(Ncells))
  
  
  for (j in 1:Ncells)
  {
    temp=sig2esti-values[j]
    counter =(abs(temp) <= Hepsi)
    histo[j]=sum(counter)/M      # Histogram Estimator
    p[j]=histo[j]/(2*Hepsi)      # PDF estimator
  }
  
  p=t(p)
  values=t(values)
  
  res<-cbind(p,values)
  return(res)
}