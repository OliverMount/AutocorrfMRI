# Plotting the axial slices as in Fig 8 of the manuscript

rm(list=ls())   # Remove all variables from workspace
library(forecast)  # For time series analysis 
library(BHMSMAfMRI) # For handling fMRI dataset
library(signal)  # for Windows and filters
library(phonTools)  # Sinc function
library(data.table)  # For shifting
library(DescTools)  # Fisher Z score
library(SuppDists)  # For Gamma distribution test for Rho
library(Matrix)  # For quadratic form
library(colorRamps)

source('TimeAC.R')
source('TimeAC0.R')
source('TimeCC.R')
source('/media/olive/Data/Rfunc/MNI2IJK.R')
source('/media/olive/Data/Rfunc/IJK2MNI.R')
source('/media/olive/Data/Rfunc/BFplot.R')
source('/media/olive/Data/Rfunc/AFplot.R')
source('/media/olive/Data/Rfunc/BF_plot.R')
source('/media/olive/Data/Rfunc/AF_plot.R')
source('/media/olive/Data/Rfunc/Fisherplot1.R')
source('/media/olive/Data/Rfunc/BFAF_plot.R')



Subpath="/mnt/NAS/fMRI/Rstate/NKIQpassed/"
flist<-list.files(Subpath)
NoS=length(flist)

Desiredq=0.05   # Desired FDR level
Zthr=5     # Z-threshold Value

# 9 mins data (Initial 5 volumes removed)
TR=c(0.645,1.4,2.5) 
Ni=c(895,399,115)      # length of time series change according to the TRs
nTR=length(TR)
d=c(60,72,60)   # Dimension of the data cube

# ROI=c("PCC","IPS","PVC","PAC","PMC","Puta","IFG")
ROI=c("PCC")
nR=length(ROI)

# ROIs coordinates (LPI orientations)
PCC=c(-2,-54,26)
IPS=c(-23,-70,46)
PVC=c(-2,-82,4)
PAC=c(-48,-24,9)
PMC=c(-38,-22,60)
Puta=c(25,-1,0)     # Putamen
IFG=c(50,23,2)

MNICoord  <- rbind(PCC,IPS,PVC,PAC,PMC,Puta,IFG)

# AFNI Filter specifications (use TRBW_MB.sh for getting the filter impulse responses)
highend=0.775
f_high=seq(0.1,highend,by=0.05)  # (0.1(1) 0.15(2) 0.2(3) 0.25(4) 0.3(5) 0.35(6) 0.4(7) 0.45(8) 0.5(9))
f_low=0.009
BW=2*(f_high-f_low)

l=1  # TR index
ns=5  # subject index
n=1  # ROI index
m=1  # BW index
N=Ni[l]


# Loading the fMRI data

fpath=paste0(Subpath,flist[ns],"/Func",l-1,"/")
BM <-drop(read.fmridata(fpath, "Afni", "10MaskedBrain+tlrc",nimages=0,dim.image=d))
NoX=sum((BM==1)*1)  # No. of voxels 


BD1 <- read.fmridata("/mnt/NAS/fMRI/templates", "Afni", "MNI152_LPI3mm+tlrc",nimages=0,dim.image=d)
BM1<-read.fmridata(fpath, "Afni", "10MaskedBrain+tlrc",nimages=0,dim.image=d)
BD=BM1*BD1
Timept=1

MNIinfo <-read.AFNI('/mnt/NAS/fMRI/templates/MNI152_LPI3mm+tlrc.HEAD')  # reading without filtering data
a=(MNIinfo$header$IJK_TO_DICOM)  # Extract the transformation matrix

#  Filter specifications  
highend=(1/(2*TR[l]))
f_high=seq(0.1,highend,by=0.05)   
f_low=0.009
BW=2*(f_high-f_low)
nBW=length(BW)

filenm=paste0(ROI[n],"_AF_",BW[m],"+tlrc.BRIK")
C_AF <- drop(read.fmridata(fpath, "Afni",filenm,nimages=0,dim.image=d))
C_BF <- drop(read.fmridata(fpath, "Afni", "PCC_BF+tlrc",nimages=0,dim.image=d))


############ Before filtering ################

# No threshold
ZNO_BF=FisherZ(C_BF)
ZNO_BF[IsZero(ZNO_BF)]=NA

# Nominal
DOFnominal=N-2
Zstat=FisherZ(C_BF)/sqrt(1/DOFnominal)
pvals=pnorm(abs(Zstat),lower.tail = F)*2 
sPvals=sort(pvals[BM==1])
Pthr=sPvals[length(sPvals[sPvals  <= Desiredq*(1:NoX)/NoX])]
Pmap=((pvals <= Pthr)*1)
Zmap_BF=(Zstat*Pmap)
Zmap_BF[IsZero(Zmap_BF)]=NA

T_Zmap_BF=((abs(Zmap_BF)>=Zthr)*1)*Zmap_BF
T_Zmap_BF[IsZero(T_Zmap_BF)]=NA


# Only Signal correlation

load(paste0("/mnt/NAS/fMRI/Rstate/NKI_DOF_BF3/",flist[ns],"_",ROI[n],"_",TR[l],".RData")) 
DOF[is.nan(DOF)]=0
DOF[DOF<=3]=0

Zstat=FisherZ(C_BF)/sqrt(1/DOF)
pvals=pnorm(abs(Zstat),lower.tail = F)*2 
sPvals=sort(pvals[BM==1])
Pthr=sPvals[length(sPvals[sPvals  <= Desiredq*(1:NoX)/NoX])]
Pmap=((pvals <= Pthr)*1)
Zmap_S=(Zstat*Pmap)
Zmap_S[IsZero(Zmap_S)]=NA

T_Zmap_S=((abs(Zmap_S)>=Zthr)*1)*Zmap_S
T_Zmap_S[IsZero(T_Zmap_S)]=NA


############## After filtering ####################

# No threshold
ZstatNO=FisherZ(C_AF)
ZstatNO[IsZero(ZstatNO)]=NA

# Nominal
DOFnominal=N-2
Zstat=FisherZ(C_AF)/sqrt(1/DOFnominal)
pvals=pnorm(abs(Zstat),lower.tail = F)*2 
sPvals=sort(pvals[BM==1])
Pthr=sPvals[length(sPvals[sPvals  <= Desiredq*(1:NoX)/NoX])]
Pmap=((pvals <= Pthr)*1)
Zmap_Nominal=(Zstat*Pmap)
Zmap_Nominal[IsZero(Zmap_Nominal)]=NA

T_Zmap_Nominal=((abs(Zmap_Nominal)>=Zthr)*1)*Zmap_Nominal
T_Zmap_Nominal[IsZero(T_Zmap_Nominal)]=NA

# Only  filter 
pathfilter=paste0("/mnt/NAS/fMRI/Rstate/NKI_filter/BP_AFNI/ImpRes_",TR[l],"_",BW[m],".1D")
g <-unlist(read.table(pathfilter,header = FALSE))
DOFfilter= (N/(sum(TimeAC0(g)^2)/(sum(abs(g)^2)^2)))

Zstat=FisherZ(C_AF)/sqrt(1/DOFfilter)
pvals=pnorm(abs(Zstat),lower.tail = F)*2 
sPvals=sort(pvals[BM==1])
Pthr=sPvals[length(sPvals[sPvals  <= Desiredq*(1:NoX)/NoX])]
Pmap=((pvals <= Pthr)*1)
Zmap_F=(Zstat*Pmap)
Zmap_F[IsZero(Zmap_F)]=NA

T_Zmap_F=((abs(Zmap_F)>=Zthr)*1)*Zmap_F
T_Zmap_F[IsZero(T_Zmap_F)]=NA


# Filter + Signal 

load(paste0("/mnt/NAS/fMRI/Rstate/NKI_DOF_AF_BU3/",flist[ns],"_",ROI[n],"_",TR[l],"_",BW[m],".RData")) 
DOF[is.nan(DOF)]=0
DOF[DOF<=3]=0

Zstat=FisherZ(C_AF)/sqrt(1/DOF)
pvals=pnorm(abs(Zstat),lower.tail = F)*2 
sPvals=sort(pvals[BM==1])
Pthr=sPvals[length(sPvals[sPvals  <= Desiredq*(1:NoX)/NoX])]
Pmap=((pvals <= Pthr)*1)
Zmap_FS=(Zstat*Pmap)
Zmap_FS[IsZero(Zmap_FS)]=NA

T_Zmap_FS=((abs(Zmap_FS)>=Zthr)*1)*Zmap_FS
T_Zmap_FS[IsZero(T_Zmap_FS)]=NA


######################     Plotting slices   ######################

sliceno=32

source('Fisherplot2.R')
Fisherplot2(BD,sliceno,Timept,ZNO_BF,ZstatNO,"Axial")   # It contains the saving option as well

source('BFAF_plot2.R')
BFAF_plot2(BD,sliceno,Timept,T_Zmap_BF,T_Zmap_S,T_Zmap_Nominal,T_Zmap_F,T_Zmap_FS,"Axial",ns)









