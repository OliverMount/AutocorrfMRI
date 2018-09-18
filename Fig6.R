# Historgram of variance of SPC in Figure 6 of the manuscript  

rm(list=ls())   # Remove all variables from workspace
library(forecast)  # For time series analysis 
library(BHMSMAfMRI) # For handling fMRI dataset
library(data.table)  # For shifting
library(DescTools)  # Fisher Z score
library(mosaic)   # Z- score
library(Matrix)  # For quadratic form
library(colorRamps)
library(devtools)
library(fields)
library(Cairo)  # For semitransperany in saving the EPS files
library(cairoDevice)
library(plyr)

# No of participants considered
flist<-c("A00058685","A00060280_2","A00060280_2R","A00060302","A00060662") 
NoS=length(flist)

TR=c(0.645,1.4,2.5) 
N=c(895,399,115)      # length of time series change according to the TRs

d=c(60,72,60)   # Dimension of the data cube (MNI voxels)
nTR=length(TR)

f_high=c(0.1,0.15,0.2,0.25,0.3)
f_low=0.009
BW=2*(f_high-f_low)

#ROI=c("PCC","IPS","PVC","PAC","PMC","Puta","IFG")
ROI=c("PCC")
nR=length(ROI)

# Choosing Subject and the ROI here
ns=4
n=1

#------------------------------------------------------------------------
#          Before Filtering histrograms (Variable TRs)
# ------------------------------------------------------------------------

l=1  # variable for TR

load(paste0("NKI_DOF_BF/",flist[ns],"_",ROI[n],"_",TR[l],".RData"))
nVAR_bf1=1/DOF
nVAR_bf1[is.infinite(nVAR_bf1)]=NA
nVAR_bf1[is.nan(nVAR_bf1)]=NA
p1<-hist(nVAR_bf1,breaks=100)
plot(p1)
dim(nVAR_bf1)<-prod(d)

l=2
load(paste0("NKI_DOF_BF/",flist[ns],"_",ROI[n],"_",TR[l],".RData"))

nVAR_bf2=1/DOF
nVAR_bf2[is.infinite(nVAR_bf2)]=NA
nVAR_bf2[is.nan(nVAR_bf2)]=NA
p2<-hist(nVAR_bf2,breaks=100)
plot(p2)
dim(nVAR_bf2)<-prod(d)

l=3
load(paste0("NKI_DOF_BF/",flist[ns],"_",ROI[n],"_",TR[l],".RData"))

nVAR_bf3=1/DOF
nVAR_bf3[is.infinite(nVAR_bf3)]=NA
nVAR_bf3[is.nan(nVAR_bf3)]=NA
p3<-hist(nVAR_bf3,breaks=100)
plot(p3)
dim(nVAR_bf3)<-prod(d)

ymax=7000
xmax=0.008
plot(p1,col="red",xlim=range(c(0,xmax)),ylim=c(0,ymax),
     ylab="Voxel counts",xlab = "Normalized variance of SPC",main="",
     cex.lab=1.5,cex.axis=1.5,cex.main=2)
plot(p2,col="green",xlim=range(c(0,1)),add=T)
plot(p3,col="blue",xlim=range(c(0,1)),add=T)
legend(3,6000, legend=c("0.645s", "1.4s","2.5s"),title="TR",
       col=c("red","green","blue"), lty=c(1,1), cex=1.2,box.lty=1)

# Data frame making
TRs=rep(TR,each=prod(d))
cVAR=c(nVAR_bf1,nVAR_bf2,nVAR_bf3)   # Combined variances
df<-data.frame(tr=TRs,variances=cVAR)
df$tr <- as.factor(df$tr)     # This is very important step 

xmin=0;xmax=0.017;xdel=0.005;
ymin=0;ymax=900;ydel=100;

par(mar=c(4,2,0,0),mfrow=c(1,1))
ggplot(df, aes(x=variances,fill=tr))+theme_bw()+ 
  geom_density(alpha=0.7)+
  scale_x_continuous(name="Variance of SPC", limits=c(xmin,xmax),breaks=seq(xmin,xmax,by=xdel)) +
  scale_y_continuous(name="Density", limits=c(ymin, ymax),breaks=seq(ymin,ymax,by=ydel))+
  scale_fill_manual(name = "TR", values = c("red", "green","blue"),labels = c("0.645 s", "1.4 s","2.5 s"))+
  theme(text = element_text(size=22),
        axis.title.x = element_text(size=20, colour="black"),
        axis.title.y = element_text(size=20, colour="black"), 
        legend.position=c(0.8,0.85),  
        plot.background = element_rect(fill="white"),
        legend.key.width = unit(0.8,"cm"),
        legend.title.align=0.5,
        legend.text.align = 0.5,
        axis.line = element_line(colour = "darkblue",size = 0.5,linetype = "solid")) +
  annotate(geom = "segment", x = (1/(N[3]-3)), xend = (1/(N[3]-3)), y = 0,
           yend = 200, color = "blue",size=1.2)+
  annotate(geom = "segment", x = (1/(N[2]-3)), xend = (1/(N[2]-3)), y = 0,
           yend = 450, color = "green",size=1.2)+
  annotate(geom = "segment", x = (1/(N[1]-3)), xend = (1/(N[1]-3)), y = 0,
           yend = 750, color = "red",size=1.2)+
  annotate(geom = "segment", x=0.0075,xend=0.01,y=600,yend=600,color = "black",size=1.2)+
  annotate("text", x = 0.014, y = 600, label = "Nominal variance",col="black",size=6)+ 
  #ggtitle("Unfiltered") +
  #theme(plot.title = element_text(hjust = 0.5))+
  ggsave(filename = "/home/olive/Desktop/TexFiles/HistBF.eps",
         device = cairo_ps,width = 6.5, height = 6)

# Mean for each density
#mu <- ddply(df, "tr", summarise, grp.mean=mean(variances,na.rm=T))
#head(mu)

# realative change from mean
#dev=c(0.001984087,0.003182739,0.008555327)
#(dev-(1/(N-3)))/dev


#------------------------------------------------------------------------
#          After Filtering histrograms
# ------------------------------------------------------------------------

################ Fixed BW, variable TRs   ############
m=1  # variable for bandwidth

l=1 # variable for TR
load(paste0("NKI_DOF_AF/",flist[ns],"_",ROI[n],"_",TR[l],"_",BW[m],".RData"))
VAR_af1=1/DOF
VAR_af1[is.infinite(VAR_af1)]=NA
VAR_af1[is.nan(VAR_af1)]=NA
p1<-hist(VAR_af1,breaks=100)
plot(p1)
dim(VAR_af1)<-prod(d)

l=2
load(paste0("NKI_DOF_AF/",flist[ns],"_",ROI[n],"_",TR[l],"_",BW[m],".RData"))
VAR_af2=1/DOF
VAR_af2[is.infinite(VAR_af2)]=NA
VAR_af2[is.nan(VAR_af2)]=NA
p2<-hist(VAR_af2,breaks=100)
plot(p2)
dim(VAR_af2)<-prod(d)

l=3
load(paste0("NKI_DOF_AF/",flist[ns],"_",ROI[n],"_",TR[l],"_",BW[m],".RData"))
VAR_af3=1/DOF
VAR_af3[is.infinite(VAR_af3)]=NA
VAR_af3[is.nan(VAR_af3)]=NA
p3<-hist(VAR_af3,breaks=100)
plot(p3)
dim(VAR_af3)<-prod(d)


############ Fixed BW and variable TRs ##########

# Data frame making
TRs=rep(TR,each=prod(d))
cVAR=c(VAR_af1,VAR_af2,VAR_af3)   # Combined variances
df<-data.frame(tr=TRs,variances=cVAR)
df$tr <- as.factor(df$tr)     # This is very important step 

xmin=0;xmax=0.030;xdel=0.005;
ymin=0;ymax=350;ydel=100;

ggplot(df, aes(x=variances,fill=tr))+theme_bw()+ 
  geom_density(alpha=0.7)+
  scale_x_continuous(name="Variance of SPC", limits=c(xmin,xmax),breaks=seq(xmin,xmax,by=xdel)) +
  scale_y_continuous(name="Density", limits=c(ymin, ymax),breaks=seq(ymin,ymax,by=ydel))+
  scale_fill_manual(name = "TR", values = c("red", "green","blue"),labels = c("0.645 s", "1.4 s","2.5 s"))+
  theme(text = element_text(size=22),
        axis.title.x = element_text(size=20, colour="black"),
        axis.title.y = element_text(size=20, colour="black"), 
        legend.position=c(0.8,0.85),  
        plot.background = element_rect(fill="white"),
        legend.key.width = unit(0.8,"cm"),
        legend.title.align=0.5,
        legend.text.align = 0.5,
        axis.line = element_line(colour = "darkblue",size = 0.5,linetype = "solid")) +
  annotate(geom = "segment", x = (1/(N[3]-3)), xend = (1/(N[3]-3)), y = 0,
           yend = 200, color = "blue",size=1.2)+
  annotate(geom = "segment", x = (1/(N[2]-3)), xend = (1/(N[2]-3)), y = 0,
           yend = 200, color = "green",size=1.2)+
  annotate(geom = "segment", x = (1/(N[1]-3)), xend = (1/(N[1]-3)), y = 0,
           yend = 200, color = "red",size=1.2)+
  annotate(geom = "segment", x=0.016,xend=0.019,y=220,yend=220,color = "black",size=1.2)+
  annotate("text", x = 0.025, y = 220, label = "Nominal variance",col="black",size=6) +


# mu <- ddply(df, "tr", summarise, grp.mean=mean(variances,na.rm=T))
# head(mu)
# dev=c(0.008400792,0.008448191,0.013460627)
# (dev-(1/(N-2)))/dev   # Realtive increase in variance

############ Fixed TR and variable BW ##########

l=1

m=1
load(paste0("NKI_DOF_AF/",flist[ns],"_",ROI[n],"_",TR[l],"_",BW[m],".RData"))
VAR_af1=1/DOF
VAR_af1[is.infinite(VAR_af1)]=NA
VAR_af1[is.nan(VAR_af1)]=NA
p1<-hist(VAR_af1,breaks=100)
plot(p1)
dim(VAR_af1)<-prod(d)

m=2
load(paste0("NKI_DOF_AF/",flist[ns],"_",ROI[n],"_",TR[l],"_",BW[m],".RData"))
VAR_af2=1/DOF
VAR_af2[is.infinite(VAR_af2)]=NA
VAR_af2[is.nan(VAR_af2)]=NA
p2<-hist(VAR_af2,breaks=100)
plot(p2)
dim(VAR_af2)<-prod(d)

m=3
load(paste0("NKI_DOF_AF/",flist[ns],"_",ROI[n],"_",TR[l],"_",BW[m],".RData"))
VAR_af3=1/DOF
VAR_af3[is.infinite(VAR_af3)]=NA
VAR_af3[is.nan(VAR_af3)]=NA
p3<-hist(VAR_af3,breaks=100)
plot(p3)
dim(VAR_af3)<-prod(d)

m=4
load(paste0("NKI_DOF_AF/",flist[ns],"_",ROI[n],"_",TR[l],"_",BW[m],".RData"))
VAR_af4=1/DOF
VAR_af4[is.infinite(VAR_af4)]=NA
VAR_af4[is.nan(VAR_af4)]=NA
p4<-hist(VAR_af4,breaks=100)
plot(p4)
dim(VAR_af4)<-prod(d)

m=5
load(paste0("NKI_DOF_AF/",flist[ns],"_",ROI[n],"_",TR[l],"_",BW[m],".RData"))
VAR_af5=1/DOF
VAR_af5[is.infinite(VAR_af5)]=NA
VAR_af5[is.nan(VAR_af5)]=NA
p5<-hist(VAR_af5,breaks=100)
plot(p5)
dim(VAR_af5)<-prod(d)



# Data frame making
BWs=rep(seq(0.1,0.3,by=0.05),each=prod(d))
cVAR=c(VAR_af1,VAR_af2,VAR_af3,VAR_af4,VAR_af5)   # Combined variances
df<-data.frame(bw=BWs,variances=cVAR)
df$bw <- as.factor(df$bw)     # This is very important step 

xmin=0;xmax=0.030;xdel=0.005;
ymin=0;ymax=350;ydel=50;

ggplot(df, aes(x=variances,fill=bw))+theme_bw()+ 
  geom_density(alpha=0.7)+
  scale_x_continuous(name="Variance of SPC", limits=c(xmin,xmax),breaks=seq(xmin,xmax,by=xdel)) +
  scale_y_continuous(name="Density", limits=c(ymin, ymax),breaks=seq(ymin,ymax,by=ydel))+
  scale_fill_manual(name = "HCF", values = c("red", "green","blue","cyan","magenta"),
                    labels = c("0.1 Hz", "0.15 Hz","0.2 Hz","0.25 Hz","0.3 Hz"))+
  theme(text = element_text(size=22),
        axis.title.x = element_text(size=20, colour="black"),
        axis.title.y = element_text(size=20, colour="black"), 
        legend.position=c(0.8,0.8),  
        plot.background = element_rect(fill="white"),
        legend.key.width = unit(0.8,"cm"),
        legend.title.align=0.5,
        legend.text.align = 0.5,
        axis.line = element_line(colour = "darkblue",size = 0.5,linetype = "solid")) +
  annotate(geom = "segment", x = (1/(N[1]-3)), xend = (1/(N[1]-3)), y = 0,
           yend = 350, color = "red",size=1.2)+
  annotate(geom = "segment", x=0.015,xend=0.019,y=180,yend=180,color = "black",size=1.2)+
  annotate("text", x = 0.025, y = 180, label = "Nominal variance",col="black",size=6)

