# Plotting the variance of SPC due to filter alone

rm(list=ls())   # Remove all variables from workspace
library(forecast)  # For time series analysis 
library(BHMSMAfMRI) # For handling fMRI dataset
library(plotly)  # For plotting
library(signal)  # for Windows and filters
library(data.table)  # For shifting
library(ggplot2)
library(reshape2)
source('TimeAC0.R')


TR=c(0.645,1.4,2.5) 
N=c(895,399,115)      # length of time series change according to the TRs

# AFNI Filter specifications 
highend=0.775
f_high=seq(0.1,highend,by=0.05)
f_low=0.009
BW=2*(f_high-f_low)

nBW=length(BW)
nTR=length(TR)

Var_AFNI<-Var_AFNI1<-Var_AFNI2<-Var_FSL<-array(NA,dim=c(nBW,nTR))  

for (l in 1:nTR) {
  temp<-array(0,dim=c((2*N[l])-1))
  for (m in 1:nBW) {
    
    pathfilter=paste0("Filters/ImpRes_",TR[l],"_",BW[m],".1D")
    if (file.exists(pathfilter))
    {
      
      # Time domain
      g <-unlist(read.table(pathfilter,header = FALSE))
      Var_AFNI[m,l]=((sum(TimeAC0(g)^2)/(sum(abs(g)^2)^2)))/N[l]

    }  
    
  }}


###################### Data frame Construction ##################

dd<-data.frame(a= f_high, b= Var_AFNI[,1],c=Var_AFNI[,2],d=Var_AFNI[,3])
colnames(dd)<-c("BW",0.645,1.4,2.5)
dd <- melt(dd ,  id.vars = 'BW', variable.name = 'TR')


# For the traingle legends in the figure

a=1/(2*TR)
b=1/(N-2)
dp<-data.frame(a=a,b=b)

a=0.5
b= 0.0115
dp1<-data.frame(a=a,b=b)


####################### Plotting the data frame ############################


xmin=0.1; xmax=0.8; xdel=0.1
ymin=0; ymax=0.017; ydel=0.002

par(mar=c(4,4.5,0,0),mfrow=c(1,1))
ggplot(dd,aes(color=TR)) + theme_bw(base_size = 30)+
  geom_point(aes(x=BW,y=value))+
  geom_line(aes(x=BW,y=value))+
  geom_point(data=dp, aes(x=a, y=b), shape=24,colour="darkorchid", size=4,fill="darkorchid") +
  geom_point(data=dp1, aes(x=a, y=b), shape=24,colour="darkorchid", size=4,fill="darkorchid") +
  scale_x_continuous(name="High Cut-off Frequency (HCF), Hz", limits=c(xmin,xmax),breaks=seq(xmin,xmax,by=xdel)) +
  scale_y_continuous(name="Variance of SPC", limits=c(ymin, ymax),breaks=seq(ymin,ymax,by=ydel))+
  annotate("text", x = 0.67, y = 0.0115, label = "Nominal variance",col="black",size=7)+
  scale_color_manual(name = "TR", values = c("red", "green","blue"),labels = c("0.645 sec", "1.4 sec","2.5 sec"))+
  theme(text = element_text(size=22),
        axis.title.x = element_text(size=20, colour="black"),
        axis.title.y = element_text(size=20, colour="black"),
        legend.position=c(0.7,0.85),
        plot.background = element_rect(fill="white"),
        legend.key.width = unit(2,"cm"),
        legend.title.align=0.5,
        legend.text.align = 0.5,
        axis.line = element_line(colour = "darkblue",size = 0.5,linetype = "solid")) 
