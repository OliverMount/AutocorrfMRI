# Z score simulation for Fig. 3a

rm(list=ls())   # Remove all variables from workspace
library(forecast)  # For time series analysis 
library(BHMSMAfMRI) # For handling fMRI dataset
library(ggplot2)  # For plotting
library(signal)  # for Windows and filters
library(phonTools)  # Sinc function
library(data.table)  # For shifting
library(corrplot)  # For plotting correlation matrix
library(mosaic)   # Z- score
library(neuRosim)  # Package for simulating fMRI data
library(MASS)    # Multivariate ???
library(devtools)
library(Hmisc)  # for fancy Violin plot
library(DescTools) # FisherZ
library(binhf)
set.seed(4776)

source('TimeAC.R')
source('TimeAC0.R')
source('TimeCC.R')
source('wantpdf.R')
source('DOFofAR1.R')
source('VarFilt.R')


True_rho <-c(0)
L=length(True_rho)
mu=c(0,0)
NN=50000  # number of times to repeat the experiment for each rho
#-------------------Filter  Specification ------------------
#TR=c(0.645,1.4,2.5)   # For NKI data
# N=c(895,399,115)      
TR=0.645
N=895  # That corresponds to that TR
# Single filter is used
f_high=0.1
f_low=0.009
BW=2*(f_high-f_low)

l=1
m=1

# Practical filters
pathfilter=paste0("ImpRes_",TR[l],"_",BW[m],".1D")
if (file.exists(pathfilter))
{
  h <-unlist(read.table(pathfilter,header = FALSE))
}



h=shift(h,((length(h)-1)/2)+1,dir="right")  # Shifting the impulse response to the center
Ch=TimeAC0(h)    # time autocorrelation of filter impulse response

#------------------Computing the Theoretical PDFs------------------
resol=0.01
rhowz=seq(-5,5,by=resol)  # Values for which PDF has to be calculated 
StdNorm=(1/sqrt(2*pi))*exp(-(resol^2)/2)  # N(0,1) distribution

#-------------Parameters for serial correlation-------------
var1=1  # Variance of process 1
var2=1  # Variance of process 2

# AR(1) parameters
alpha=0.4
beta=0.9


Tpdf<-TpdfZ<-array(0,dim=c(length(rhowz),L))
me<-va<-array(0,dim=c(L))

# Initialization
AFcorr<-AFcorrZ<-BF_IID<-BF_IIDZ<-BFcorr<-BFcorrZ<-array(0,dim=c(NN,L))   # After filtering

for (k in 1:L)  # length rho
{
  
  cmtx=rbind(c(var1, True_rho*sqrt(var1)*sqrt(var2)),
             c(True_rho*sqrt(var1)*sqrt(var2), var2))
  for (kk in 1:NN)  # length of trails
  {
    
    r<-mvrnorm(N,mu,Sigma=cmtx)  # Bivariate Gaussian (No serial correlation)
    
    # Correlation Before filtering (No AR(1))
    BF_IID[kk,k]=cor(r[,1],r[,2])  # raw
    BF_IIDZ[kk,k]=FisherZ(BF_IID[kk,k])  # Fisher Z
    
    
    if (alpha !=0 & beta !=0)
    {
      #--------- Generate two AR(1) processes -----------------
      v <- array(0, dim = c(N)) # Initializing
      w <- array(0, dim = c(N)) # Initializing
      
      # Serially correlated sequence
      for (i in 1:(N-1))
      {
        v[i+1]=(alpha*v[i])+r[i+1,1]  # x and y are AR(1) Processes
        w[i+1]=(beta*w[i])+r[i+1,2]
      }
    }
    else
    {
      v=r[,1]
      w=r[,2]
    }
    
    # Correlation Before filtering (With AR(1))
    BFcorr[kk,k]=cor(v,w)  # raw
    BFcorrZ[kk,k]=FisherZ(BFcorr[kk,k])  # Zscore
    
    
    #------------ Bandpass filtering ----------------------
    
    x<-convolve(v[1:N],rev(h),type="o")   # should be linear convolution 
    y<-convolve(w[1:N],rev(h),type="o")  
    
    # Correlation after filtering
    AFcorr[kk,k]=cor(x,y)  # raw
    AFcorrZ[kk,k]=FisherZ(AFcorr[kk,k])  # Zscore
  }
}

#------------Computing the Simulation PDFs (After filtering -------------

binsize=40
IID<-NC_BF <-NC_AF <-C_BF <-C_F <- C_FS <-array(0,dim=c(binsize,binsize))

########### No Correction part ############

temp=VarFilt(Ch,alpha,beta,N,True_rho,var1,var2) # it returns mean and variance
va_filter=((sum(Ch)^2)/((sum(abs(h)^2))^2))/N

Ch1=TimeAC0(c(rep(0,(N-1)/2),1,rep(0,(N-1)/2)))
temp1=DOFofAR1(Ch1,alpha,beta,N,True_rho,var1,var2)  # it returns mean and variance

# IID
Noco1=(BF_IIDZ-FisherZ(True_rho))/sqrt(1/(N-3))
IID=wantpdf(Noco1,binsize)

# BF  (Nominal)
Noco2=(BFcorrZ-FisherZ(temp1[1]))/sqrt(1/(N-3))
NC_BF=wantpdf(Noco2,binsize)

# AF (Nominal)
Noco3=(AFcorrZ-FisherZ(temp[1]))/sqrt(1/(N-3))
NC_AF=wantpdf(Noco3,binsize)

######### Corrections part ###############
# BF (Signal only correction)
sc=(1/(1-(temp1[1]^2)))^2
Correction0=(BFcorrZ-FisherZ(temp1[1]))/sqrt(temp1[2]*sc)
C_BF=wantpdf(Correction0,binsize) 

# Only Filter correction
Correction1=(AFcorrZ-FisherZ(temp[1]))/sqrt(va_filter)
C_F=wantpdf(Correction1,binsize)

# Filter + Signal correction
sc=(1/(1-(temp[1]^2)))^2
Correction2=(AFcorrZ-FisherZ(temp[1]))/sqrt(temp[2]*sc)
C_FS=wantpdf(Correction2,binsize) 


##########---------PLOTTING THE RESULTS -----------#####

plot(NC_BF[(binsize+1):(2*binsize)],NC_BF[1:binsize],type="l",
     xlim =c(-8,8),ylim=c(0,0.7),xlab="Z-score",ylab="Probability density",cex.lab=1.5,
     xaxt="n",yaxt="n",xaxs="i", yaxs="i",lwd=2,col="blue")
lines(C_BF[(binsize+1):(2*binsize)],C_BF[1:binsize],type="l",lwd=2,lty=2,col="black")
lines(NC_AF[(binsize+1):(2*binsize)],NC_AF[1:binsize],type="l",lwd=2,col="orange")
lines(C_F[(binsize+1):(2*binsize)],C_F[1:binsize],type="l",lwd=2,col="green")
lines(C_FS[(binsize+1):(2*binsize)],C_FS[1:binsize],type="l",lwd=2,col="red")
xtick<-seq(-8,8, by=1)
axis(side=1, at=xtick,cex.axis=1.5)
ytick<-seq(0,0.7,by=0.1)
axis(side=2, at=ytick,cex.axis=1.5)
grid(NULL, NULL, col = "gray", lty = "dotted")
legend(-6, 0.7, legend=c("Unfiltered (Nominal)","Unfiltered (Signal only correction)",
                         "Filtered (Nominal)","Filtered (Filter only correction)",
                         "Filtered (Filter + Signal correction)"),
       col=c("blue","black","orange","green","red"), lty=c(1,2,1,1,1),cex=1.3,box.lty=0,lwd=c(3,3,3,3))
lines(c(1.96,1.96),c(0,.45),type="l",lwd=2,lty=2,col="black")
lines(c(-1.96,-1.96),c(0,.45),type="l",lwd=2,lty=2,col="black")
text(4,0.3,labels="Threshold",cex=1.3)




