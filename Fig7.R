# This program for calculation of FPR 

rm(list=ls())   # Remove all variables from workspace
library(forecast)  # For time series analysis 
library(BHMSMAfMRI) # For handling fMRI dataset
library(signal)  # for Windows and filters
library(data.table)  # For shifting
library(DescTools)  # Fisher Z score
library(SuppDists)  # For Gamma distribution test for Rho
library(Matrix)  # For quadratic form
library(colorRamps)
library(devtools)
library(plyr)
library(latex2exp)

source('/media/olive/Data/Rfunc/TimeAC0.R')
source('/media/olive/Data/Rfunc/TimeAC.R')
source('/media/olive/Data/Rfunc/TimeCC.R')
source('/media/olive/Data/Rfunc/VoxelDOF.R')

# get the subjects list 

flist<-c("A00058685","A00060280_2","A00060280_2R","A00060302","A00060662") 
NoS=length(flist)

# Choose the TR and the BW to plot 
# 9 mins data (Initial 5 volumes removed)
TR=c(0.645,1.4,2.5) 
Ni=c(895,399,115)      # length of time series change according to the TRs

d=c(60,72,60)   # Dimension of the data cube

nTR=length(TR)

# ROI=c("PCC","IPS","PVC","PAC","PMC","Puta","IFG")
ROI=c("PCC")
nR=length(ROI)

# ROIs coordinates (LPI coordinate order)
PCC=c(-2,-54,26)

tri  <- matrix(c(PCC,IPS,PVC,PAC,PMC,Puta,IFG),nrow=length(ROI))
NAF_F<-NAF_FS<-array(0,dim=c(20,nTR,nR,NoS))
NAF_conv<-array(0,dim=c(nR,NoS))

nTR=length(TR)

FPR_BF<-FPR_BF_new<-array(NA,dim=c(nTR,nR,NoS)) # Storing BF FPR
FPR_AF<-FPR_AF_new0<-FPR_AF_new1<-FPR_AF_new2<-array(NA,dim=c(14,nTR,nR,NoS))  # Storing AF FPR
dof<-array(NA,dim=c(14,nTR)) 

# desired q-value
Desiredq=0.001

for (ns in 1:NoS) {
  for (l in 1:nTR) {
    
    #------------------------------------------------------------------------
    #        False positive rate (FPR)  Before Filtering  
    # ------------------------------------------------------------------------
    
    # Loading required parameters
    
    
    fpath=paste0(flist[ns],"/Func",l-1,"/")
    BM <-drop(read.fmridata(fpath, "Afni", "10MaskedBrain+tlrc",nimages=0,dim.image=d))
    NoX=sum((BM==1)*1)  # No. of voxels
    
    
    #  Filter specifications  
    highend=(1/(2*TR[l]))
    f_high=seq(0.1,highend,by=0.05)   
    f_low=0.009
    BW=2*(f_high-f_low)
    nBW=length(BW)
    
    N=Ni[l]
    
    #------------------------------------------------------------------------
    #             False positive rate (FPR)     Before Filtering 
    # ------------------------------------------------------------------------
    
    for (m in 1:nBW) {
      
      for (n in 1:nR)  # No. of RSN
      {
        filenm=paste0(ROI[n],"_BF","+tlrc.BRIK")
        C_BF <- drop(read.fmridata(fpath, "Afni", filenm,nimages=0,dim.image=d))
        
        ##   Nominal
        Zstat=FisherZ(C_BF)/sqrt(1/(N-2))
        pvals=pnorm(abs(Zstat),lower.tail = F)*2
        sPvals=sort(pvals[BM==1])
        Pthr=sPvals[length(sPvals[sPvals  <= Desiredq*(1:NoX)/NoX])]
        FPR_BF[l,n,ns]<- (sum(((pvals<=Pthr)*1))/NoX)
        
        ##   Signal correlation
        
        load(paste0("NKI_DOF_BF/",flist[ns],"_",ROI[n],"_",TR[l],".RData"))
        DOF[is.nan(DOF)]=0
        DOF[DOF<=3]=0
        
        Zstat=FisherZ(C_BF)/sqrt(1/DOF)
        pvals=pnorm(abs(Zstat),lower.tail = F)*2
        sPvals=sort(pvals[BM==1])
        Pthr=sPvals[length(sPvals[sPvals  <= Desiredq*(1:NoX)/NoX])]
        FPR_BF_new[l,n,ns]<- (sum(((pvals<=Pthr)*1))/NoX)
        
        #------------------------------------------------------------------------
        #             False positive rate (FPR)     After Filtering 
        # ------------------------------------------------------------------------
        
        pathfilter=paste0("Filter/ImpRes_",TR[l],"_",BW[m],".1D")
        if (file.exists(pathfilter))
        {
          
          filenm=paste0(ROI[n],"_AF_",BW[m],"+tlrc.BRIK")
          C_AF <- drop(read.fmridata(fpath, "Afni", filenm,nimages=0,dim.image=d))
          
          
          # Nominal  (N-2)
          
          DOFnominal=N-2
          Zstat=FisherZ(C_AF)/sqrt(1/DOFnominal)
          pvals=pnorm(abs(Zstat),lower.tail = F)*2 
          sPvals=sort(pvals[BM==1])
          Pthr=sPvals[length(sPvals[sPvals  <= Desiredq*(1:NoX)/NoX])]
          FPR_AF[m,l,n,ns]=(sum(((pvals<=Pthr)*1))/NoX)
          
          # Only Filter (Method from literature)
          g <-unlist(read.table(pathfilter,header = FALSE))
          sp=fft(g)  # Filter frequency response 
          DOFfilter=(sum(abs(sp)^2)^2)/(sum(abs(sp)^4))
          
          Zstat=FisherZ(C_AF)/sqrt(1/DOFfilter)
          pvals=pnorm(abs(Zstat),lower.tail = F)*2 
          sPvals=sort(pvals[BM==1])
          Pthr=sPvals[length(sPvals[sPvals  <= Desiredq*(1:NoX)/NoX])]
          FPR_AF_new0[m,l,n,ns]=(sum(((pvals<=Pthr)*1))/NoX)
          
          
          # Only filter (our method)
          
          dof[m,l]=(N/(sum(TimeAC0(g)^2)/(sum(abs(g)^2)^2)))
          
          Zstat=FisherZ(C_AF)/sqrt(1/dof[m,l])
          pvals=pnorm(abs(Zstat),lower.tail = F)*2 
          sPvals=sort(pvals[BM==1])
          Pthr=sPvals[length(sPvals[sPvals  <= Desiredq*(1:NoX)/NoX])]
          FPR_AF_new1[m,l,n,ns]=(sum(((pvals<=Pthr)*1))/NoX)
          
          # Filter + Signal
          
          load(paste0("NKI_DOF_AF/",flist[ns],"_",ROI[n],"_",TR[l],"_",BW[m],".RData")) 
          
          DOF[is.nan(DOF)]=0
          DOF[DOF<=3]=0
          
          Zstat=FisherZ(C_AF)/sqrt(1/DOF)
          pvals=pnorm(abs(Zstat),lower.tail = F)*2 
          sPvals=sort(pvals[BM==1])
          Pthr=sPvals[length(sPvals[sPvals  <= Desiredq*(1:NoX)/NoX])]
          FPR_AF_new2[m,l,n,ns]<- (sum(((pvals<=Pthr)*1))/NoX)
        }
      }
    }}
}

########  Data frame making ###############
ROIin=1

# AFter Filtering

# Nominal
highend=0.775
a= seq(0.1,highend,by=0.05)   # Note that this starts with 0.1 (not 0.05)
bN=apply(FPR_AF[,,ROIin,],c(1,2),mean)
cN=apply(FPR_AF[,,ROIin,],c(1,2),sd)
df<-data.frame(a,bN)
colnames(df)<-c("BW",0.645,1.4,2.5)
df <- melt(df,id.vars="BW",variable.name="TR")
df$SD<-as.vector(cN/sqrt(NoS))
df$Method<-rep("Nominal",length(a))

# From Literature
bL=apply(FPR_AF_new0[,,ROIin,],c(1,2),mean)
cL=apply(FPR_AF_new0[,,ROIin,],c(1,2),sd)
df0<-data.frame(a,bL)
colnames(df0)<-c("BW",0.645,1.4,2.5)
df0 <- melt(df0,id.vars="BW",variable.name="TR")
df0$SD<-as.vector(cL/sqrt(NoS))
df0$Method<-rep("Literature",length(a))

# Only filter
bF=apply(FPR_AF_new1[,,ROIin,],c(1,2),mean)
cF=apply(FPR_AF_new1[,,ROIin,],c(1,2),sd)
df1<-data.frame(a,bF)
colnames(df1)<-c("BW",0.645,1.4,2.5)
df1 <- melt(df1,id.vars="BW",variable.name="TR")
df1$SD<-as.vector(cF/sqrt(NoS))
df1$Method<-rep("Filter only",length(a))

# filter + signal
bFS=apply(FPR_AF_new2[,,ROIin,],c(1,2),mean)
cFS=apply(FPR_AF_new2[,,ROIin,],c(1,2),sd)
df2<-data.frame(a,bFS)
colnames(df2)<-c("BW",0.645,1.4,2.5)
df2 <- melt(df2,id.vars="BW",variable.name="TR")
df2$SD<-as.vector(cFS/sqrt(NoS))
df2$Method<-rep("Filter+Signal",length(a))

##################  Nominal and Only filter #####################

ls=1.2  # Line size
Method=c("Nominal","Filter only")
ggplot(df, aes(BW,value,color=TR,linetype=Method)) + theme_bw()+
  geom_line(linetype=2,size=ls) +
  geom_point(size=ls) +
  geom_errorbar(aes(ymin=value-SD, ymax=value+SD), width=.01)+
  scale_x_continuous(name="High Cut-off Frequency (HCF), Hz", limits=c(0.1,highend),breaks=a) +
  scale_y_continuous(name="False Positive Rate", limits=c(0,0.7),breaks=seq(0,0.7,by= 0.1)) + theme_classic()+ theme(text = element_text(size=15))+
  geom_line(data=df1,linetype=1,size=ls)+ geom_point(data=df1,size=ls)   +
  geom_errorbar(data=df1,aes(ymin=value-SD, ymax=value+SD), width=.01)+
  theme(legend.position ="top" , 
        legend.justification = c("center", "top"),
        legend.box.just = "left",
        plot.background = element_rect(fill="white"), 
        legend.margin = margin(1,1,-4,1,"cm"),
        legend.title = element_text(size=20, color = "blue"),
        legend.text= element_text(size=20),
        legend.title.align=0.5,
        legend.text.align=0.5,
        legend.key.width = unit(2,"cm"),  # width of the line
        axis.title.x = element_text(size=16, colour="black"),
        axis.title.y = element_text(size=18, colour="black"),
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=15),
        text = element_text(size=15),
        panel.grid.major = element_line(size = 0.3, linetype = 'dotted',
                                        colour = "gray88"))+ 
  annotate("text", x = 0.6, y = 0.51, label = "FDR corrected, q < 0.001",col="black",cex=5)+
  labs(x = "High Cut-off Frequency (HCF), Hz", y = "False Positive Rate", linetype = "Method",color="TR")  +
  scale_linetype_manual(labels = c("Nominal","Filter only"), values = c(2,1) )+
  scale_color_manual(name = "TR", values = c("red", "green","blue"),labels = c("0.645 s", "1.4 s","2.5 s"))+
  guides(color = guide_legend(title.position = "top",order = 1, nrow = 3,
                              override.aes = list(size=6,shape=NA), keywidth=1),
         linetype = guide_legend(title.position = "top",nrow = 2, byrow = TRUE,
                                 override.aes = list(size=1)))


#######  Nominal and filter + signal ###############


ls=1.2  # Line size
ggplot(df, aes(BW,value,color=TR,linetype=Method)) +theme_bw()+ 
  geom_line(linetype=2,size=ls) +
  geom_point(size=ls) +
  geom_errorbar(aes(ymin=value-SD, ymax=value+SD), width=.01)+
  scale_x_continuous(name="High cut-off frequency (HCF), Hz", limits=c(0.1,highend),breaks=a) +
  scale_y_continuous(name="False Positive Rate", limits=c(0,0.7),breaks=seq(0,0.7,by= 0.1)) + theme_classic()+ theme(text = element_text(size=15))+
  geom_line(data=df2,linetype=1,size=ls)+ geom_point(data=df2,size=ls)   +
  geom_errorbar(data=df2,aes(ymin=value-SD, ymax=value+SD), width=.01)+
  theme(legend.position ="top" , 
        legend.justification = c("center", "top"),
        legend.box.just = "left",
        legend.background = element_rect(fill="white"),
        legend.margin = margin(1,1,-4,1,"cm"),
        legend.title = element_text(size=20, color = "blue"),
        legend.text= element_text(size=20),
        legend.title.align=0.5,
        legend.text.align = 0.5,
        legend.key.width = unit(2,"cm"),  # width of the line
        axis.title.x = element_text(size=16, colour="black"),
        axis.title.y = element_text(size=18, colour="black"),
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=15),
        text = element_text(size=15),
        panel.grid.major = element_line(size = 0.3, linetype = 'dotted',
                                        colour = "gray88"))+ 
  annotate("text", x = 0.6, y = 0.51, label = "FDR corrected, q < 0.001",col="black",cex=5)+
  labs(x = "High Cut-off Frequency (HCF), Hz", y = "False Positive Rate", linetype = "Method",color="TR")  +
  scale_linetype_manual(labels = c("Nominal","Filter + Signal "), values = c(2,1) )+
  scale_color_manual(name = "TR", values = c("red", "green","blue"),labels = c("0.645 s", "1.4 s","2.5 s"))+
  guides(color = guide_legend(title.position = "top",order = 1, nrow = 3,
                              override.aes = list(size=6,shape=NA), keywidth=1),
         linetype = guide_legend(title.position = "top",nrow = 2, byrow = TRUE,
                                 override.aes = list(size=1)))



### Data frame for only a particular frequency for making bar chart


#After Filtering
BWin=1
Methods=c("Conventional","Filter only","Filter+Signal")
a=rep(TR,length(Methods))
b=rep(Methods,each=length(TR))
c=c(bN[BWin,],bF[BWin,],bFS[BWin,])
d=c(cN[BWin,],cF[BWin,],cFS[BWin,])
bardf<-data.frame(a,b,c,d)
colnames(bardf)<-c("TR","Methods","Mean","SD")
bardf$TR=as.factor(bardf$TR)

#######  After filtering bar chart  ###############

ggplot(bardf, aes(x=TR,y=Mean, fill=Methods)) + 
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean+SD), width=.2,
                position=position_dodge(.9))+ 
  theme_classic()+
  scale_y_continuous(limits=c(0,0.6),breaks=seq(0,1,by= 0.1))+
  theme(legend.position = c(.52, 1), 
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(1,2,0,1,"cm"),
        legend.title.align=0.5,
        legend.text.align = 0.5,
        legend.title = element_text(size=20, color = "blue"),
        legend.text=element_text(size=20), 
        axis.title.x = element_text(size=16, colour="black"),
        axis.title.y = element_text(size=18, colour="black"),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=16),
        text = element_text(size=15),
        panel.grid.major = element_line(size = 0.3, linetype = 'dotted',
                                        colour = "gray88")) +
  labs(x = "TR (sec)", y = "False positive rate", fill = "Method") +
  scale_fill_manual(labels =  c("Nominal","Filter only","Filter + Signal"), values = c("red", "green","blue") )+
  annotate("text", x = 3, y = 0.35, label = "FDR corrected, q < 0.001",col="black",cex=5)

# Before Filtering

Methods=c("Nominal","Signal only")
a=rep(TR,length(Methods))
b=rep(Methods,each=length(TR))
bFN=apply(FPR_BF[,ROIin,],1,mean)# Nominal
cFN=apply(FPR_BF[,ROIin,],1,sd)
bFS=apply(FPR_BF_new[,ROIin,],1,mean)# Exploiting Signal only
cFS=apply(FPR_BF_new[,ROIin,],1,sd)
c=c(bFN,bFS)
d=c(cFN,cFS)
barBF<-data.frame(a,b,c,d)
colnames(barBF)<-c("TR","Methods","Mean","SD")
barBF$TR=as.factor(barBF$TR)


#######  Before filtering bar chart  ###############

ggplot(barBF, aes(x=TR,y=Mean, fill=Methods)) + 
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean+SD), width=.2,
                position=position_dodge(.9))+ 
  theme_classic()+
  scale_y_continuous(limits=c(0,0.6),breaks=seq(0,1,by= 0.1))+
  theme(legend.position = c(.6, 1), 
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(1,2,0,1,"cm"),
        legend.title.align=0.5,
        legend.text.align = 0.5,
        legend.title = element_text(size=20, color = "blue"),
        legend.text=element_text(size=20), 
        axis.title.x = element_text(size=16, colour="black"),
        axis.title.y = element_text(size=18, colour="black"),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=16),
        text = element_text(size=15),
        panel.grid.major = element_line(size = 0.3, linetype = 'dotted',
                                        colour = "gray88")) +
  labs(x = "TR (sec)", y = "False positive rate", fill = "Method",Main="Before Filtering") +
  scale_fill_manual(labels =  c("Nominal","Signal only"), values = c("yellow","darkorchid") )+
  annotate("text", x = 3, y = 0.41, label = "FDR corrected, q < 0.001",col="black",cex=5) 


#######  Supplementary S2(b) figure:   Filter (literaure)  and Filter (Our Method)  ###############

ggplot(df1, aes(BW,value,color=TR,linetype=Method)) +theme_bw()+ 
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=value-SD, ymax=value+SD), width=.01)+
  scale_x_continuous(name="High cut-off frequency (HCF), Hz", limits=c(0.1,highend),breaks=a) +
  scale_y_continuous(name="False Positive Rate", limits=c(0,0.5),breaks=seq(0,0.5,by= 0.1)) + theme_classic()+ theme(text = element_text(size=15))+
  geom_line(data=df0)+ geom_point(data=df0)   +
  geom_errorbar(data=df0,aes(ymin=value-SD, ymax=value+SD), width=.01)+
  theme(legend.position ="top" , 
        legend.justification = c("center", "top"),
        legend.box.just = "left",
        legend.background = element_rect(fill="white"),
        legend.margin = margin(1,1,-4,1,"cm"),
        legend.title = element_text(size=20, color = "blue"),
        legend.text= element_text(size=20),
        legend.title.align=0.5,
        legend.text.align = 0.5,
        legend.key.width = unit(2,"cm"),  # width of the line
        axis.title.x = element_text(size=20, colour="black"),
        axis.title.y = element_text(size=20, colour="black"),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        text = element_text(size=15),
        panel.grid.major = element_line(size = 0.3, linetype = 'dotted',
                                        colour = "gray"))+ 
  annotate("text", x = 0.4, y = 0.35, label = "FDR corrected, q < 0.001",col="black",cex=6)+
  labs(x = "High Cut-off Frequency (HCF), Hz", y = "False Positive Rate", linetype = "Methods",color="TR")  +
  scale_linetype_manual(labels = c("Our approach","Davey et. al [Eqn. (8)]"), values = c(1,2) )+
  scale_color_manual(name = "TR", values = c("red", "green","blue"),labels = c("0.645 s", "1.4 s","2.5 s"))+
  guides(color = guide_legend(title.position = "top",order = 1, nrow = 3,
                              override.aes = list(size=6,shape=NA), keywidth=1),
         linetype = guide_legend(title.position = "top",nrow = 2, byrow = TRUE,
                                 override.aes = list(size=1)))