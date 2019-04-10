# Simulation of AR(1) network of 10 nodes in Fig 2 and Fig 4

rm(list=ls())   # Remove all variables from workspace
library(forecast)  # For time series analysis 
library(BHMSMAfMRI) # For handling fMRI dataset
library(signal)  # for Windows and filters
library(data.table)  # For shifting
library(DescTools)  # Fisher Z score
library(Matrix)  # For quadratic form
library(colorRamps)
library(igraph)
library(qgraph)
library(rgl)
library(fields)  # image.plot
library(binhf)  # For shift function
library(s2dverification)

source('TimeAC.R')
source('TimeAC0.R')  # Not normalized
source('TimeCC.R')
source('gset.R')

Desiredq=0.001    # Desired q value (corrected p-value)
Nn=9  # Number of nodes
PstatTrue=c(-0.21,0.58,0.72,0.4,-0.1,0.8,-0.68,0.32,0.63)  # raw Pearson Statistics (True)
PstatBF=c(-0.004,0.51,0.05,0.4,-0.1,0.46,-0.645,0.29,0.63) #  raw Pearson Statistics (Before Filtering, AR correlated)
PstatAF=c(-0.0175,0.57,0.35, 0.4,-0.1, 0.48,-0.657,0.3,0.63) #  raw Pearson Statistics (After Filtering, AR and filter correlated)

vaBF=c(5.470601e-05,0.002801171,0.0003199577,0.1254442,0.02117411,0.2641356,0.0008700172,0.002305576,0.1016729) #Variance before Filtering

## Filter 
TR=1.4
N=399
f_high=0.1
f_low=0.009
BW=2*(f_high-f_low)

# Realistic AFNI filters
pathfilter=paste0("Filter/ImpRes_",TR,"_",BW,".1D")
if (file.exists(pathfilter))
{
  h <-unlist(read.table(pathfilter,header = FALSE))
}
h=binhf::shift(h,((length(h)-1)/2)+2,dir="right")
Ch=TimeAC0(h)

# GRAPH layout
coord=array(0,dim=c(10,2))
coord[1,]=c(3,3)       # Seed

coord[2,]=c(5.5,2.5)    # Group1
coord[3,]=c(6,3)
coord[4,]=c(5.5,3.5)

coord[5,]=c(3,6)      # Group2
coord[6,]=c(4,5.5)
coord[7,]=c(2,5.5)

coord[8,]=c(0,3)    # Group 3
coord[9,]=c(0.5,2.5)
coord[10,]=c(0.5,3.5)

######### Making graph in Fig.2  ########

g <- graph(edges=c(1,2,1,3,1,4,1,5,1,6,1,7,1,8,1,9,1,10),n=10,directed=F)
E(g)$weight=PstatAF
g <- gset(g)
V(g)$label.cex=1.5

zr<-range(PstatAF)
Ocolorraw=blue2red(Nn)
OcolorFF<-Ocolorraw
Ocolorraw[order(PstatAF)]<-OcolorFF
Ocolorraw<-OcolorFF[c(3,8,5,6,2,7,1,4,9)]

bre=seq(min(PstatAF),max(PstatAF),by=(max(PstatAF)-min(PstatAF))/(length(Ocolorraw)))
bre1=sort(c(round(PstatAF,2),0.65))

##### Plotting Fig 2 #######

par(mar=c(8,0,0,0),mfrow=c(1,1))
V(g)$color=c("#FF0000",Ocolorraw)
plot.igraph(g,layout=coord,mark.groups = list(c(2,3,4),c(5,6,7),c(8,9,10)),
            mark.col ="lightblue",vertex.label.family="serif",mark.border=F)
text(0,-.93,"Seed",cex = 1.5,col="blue",font=2)
#
image.plot(legend.only=TRUE,
           zlim= zr,horizontal = T,col=OcolorFF,
           legend.shrink = 0.8, legend.width = 1,
           #breaks=bre,
           lab.breaks = bre1, #c("a","b","c","d","e","f","g","h","i","j"),
           legend.lab = "Pearson Correlation",legend.cex = 1.3,
           smallplot= c(0.1,0.9,0.15,0.2))

######## For Fig. 4a and 4d ################ 
#### Preparing the graphs (plotting done afterwards)

# Before filtering (FisherZ graph)
g <- graph(edges=c(1,2,1,3,1,4,1,5,1,6,1,7,1,8,1,9,1,10),n=10,directed=F)
CovalBF=FisherZ(PstatBF)  # FsiherZ(Correlation value)
g <- gset(g)
V(g)$label.cex=3

# After filtering (FisherZ graph)
g1 <- graph(edges=c(1,2,1,3,1,4,1,5,1,6,1,7,1,8,1,9,1,10),n=10,directed=F)
CovalAF=FisherZ(PstatAF)  # FisherZ(Correlation value)
g1 <- gset(g1)
V(g1)$label.cex=3

# Dealing with colors for the combined plot
FZarr=c(CovalBF,CovalAF)
OcolorF=blue2red(length(FZarr))
OcolorFF<-OcolorF1<-OcolorF

OcolorF[order(FZarr)]<-OcolorFF   # Order(zarr) is node values
Cval=c("#FF0000",OcolorF[1:9],"#FF0000",OcolorF[10:18])  # color value
V(g)$color<-c("#FF0000",OcolorF[1:9])
V(g1)$color<-c("#FF0000",OcolorF[10:18])
FZval=sort(FZarr)  # Sorted Fisher value

### Plotting the graphs
###### Fig 4a ##########
par(mar=c(0,0,0,0),mfrow=c(1,1))
plot.igraph(g,layout=coord,mark.groups = list(c(2,3,4),c(5,6,7),c(8,9,10)),
            mark.col ="lightblue",vertex.label.family="serif",mark.border=F)
text(0,-.93,"Seed",cex = 3,col="blue",font=2)

dev.off()

###### Fig 4d ##########

par(mar=c(0,0,0,0),mfrow=c(1,1))
plot.igraph(g1,layout=coord,mark.groups = list(c(2,3,4),c(5,6,7),c(8,9,10)),
            mark.col ="lightblue",vertex.label.family="serif",mark.border=F)
text(0,-.93,"Seed",cex = 3,col="blue",font=2)

###### Color bar for Fig 4a and 4d ##########

ma=max(FZarr)
mi=min(FZarr)

scale = (length(OcolorF1))/(ma-mi)
ticks=round(seq(mi,ma, len=11),1)

plot(1,1, type='n', bty='n',xaxt='n', yaxt='n', ylab='',xlab='',xlim=c(mi,ma),ylim=c(0,5))
axis(side=1,at=ticks,cex.axis=2.5)
mtext(side=1, text="Fisher Z", line=4,cex=3)
for (i in 1:(length(OcolorF1))) {
  y = (i-1)/scale + mi
  rect(y,-1,y+1,.4,col=OcolorF1[i], border=NA)    # Change the last value here for the height of the color bar
}

######## Preparing for graphs Fig. 4b,c,e,f,g ################ 

######## Setting up the Z-value graphs (Before Filtering)
# Nominal
DOFnominal=N-3
ZstatN=FisherZ(PstatBF)/sqrt(1/DOFnominal)
pvalsN=pnorm(abs(ZstatN),lower.tail = F)*2 
Lp=length(pvalsN)
qvals=p.adjust(pvalsN, method = c("BH"), n = Lp)

sPvals=sort(pvalsN)             # Finding P-value threshold and thresholding with p-values
Pthr=sPvals[length(sPvals[sPvals <= Desiredq*(1:Nn)/Nn])]
pvalsN <= Pthr
PattN_BF=(pvalsN <= Pthr)*1
sum(((pvalsN<=Pthr)*1))


# Signal only

ZstatS=FisherZ(PstatBF)/sqrt(vaBF)
pvals=pnorm(abs(ZstatS),lower.tail = F)*2 
Lp=length(pvals)
qvals=p.adjust(pvals, method = c("BH"), n = Lp)

sPvals=sort(pvals)
Pthr=sPvals[length(sPvals[sPvals <= Desiredq*(1:Nn)/Nn])]
pvals <= Pthr
PattS=(pvals <= Pthr)*1
sum(((pvals<=Pthr)*1))

EdgN=PattN_BF*c(2:(Nn+1))
EdgN=EdgN[EdgN!=0]

EdgS=PattS*c(2:(Nn+1))
EdgS=EdgS[EdgS!=0]


######## Setting up the Z-value graphs (After Filtering case) ##############
# Nominal
DOFnominal=N-3
ZstatN_AF=FisherZ(PstatAF)/sqrt(1/DOFnominal)
pvalsN=pnorm(abs(ZstatN_AF),lower.tail = F)*2 
Lp=length(pvalsN)
qvals=p.adjust(pvalsN, method = c("BH"), n = Lp)


sPvals=sort(pvalsN)                 # Finding P-value threshold and thresholding with p-values
Pthr=sPvals[length(sPvals[sPvals <= Desiredq*(1:Nn)/Nn])]
pvalsN <= Pthr
PattN_AF=(pvalsN <= Pthr)*1
sum(((pvalsN<=Pthr)*1))


# DOF filter
DOFfilter=1/((sum(Ch^2)/((sum(abs(h)^2))^2))/N)
ZstatF=FisherZ(PstatAF)/sqrt(1/DOFfilter)
pvals=pnorm(abs(ZstatF),lower.tail = F)*2 
Lp=length(pvals)
qvals=p.adjust(pvals, method = c("BH"), n = Lp)


sPvals=sort(pvals)
Pthr=sPvals[length(sPvals[sPvals <= Desiredq*(1:Nn)/Nn])]
pvals <= Pthr
PattF=(pvals <= Pthr)*1
sum(((pvals<=Pthr)*1))

# Variable DOF
DOFeffective=c(373,110,146,38,56,8,323,130,12.5)
ZstatFS=FisherZ(PstatAF)/sqrt(1/DOFeffective)
pvals=pnorm(abs(ZstatFS),lower.tail = F)*2 
Lp=length(pvals)
qvals=p.adjust(pvals, method = c("BH"), n = Lp)

sPvals=sort(pvals)
Pthr=sPvals[length(sPvals[sPvals <= Desiredq*(1:Nn)/Nn])]
pvals <= Pthr
PattFS=(pvals <= Pthr)*1
sum(((pvals<=Pthr)*1))

EdgN_AF=PattN_AF*c(2:(Nn+1))
EdgN_AF=EdgN_AF[EdgN_AF!=0]

EdgF=PattF*c(2:(Nn+1))
EdgF=EdgF[EdgF!=0]

EdgFS=PattFS*c(2:(Nn+1))
EdgFS=EdgFS[EdgFS!=0]


# Dealing with color bars for the combined plot

Zarr=c(ZstatN,ZstatS,ZstatN_AF,ZstatF,ZstatFS)
ZL=length(Zarr)
OcolorZ=blue2red(ZL)
OcolorZZ=OcolorZ
OcolorZ1=OcolorZ

OcolorZ[order(Zarr)]<-OcolorZZ   # Order(zarr) is node values
Zval=sort(Zarr)  # Correlation value


# Nominal (Before filtering )
a=rep(1,length(EdgN))
gN <- graph(edges=c(rbind(a,EdgN)),n=10,directed=F)
gN <- gset(gN)
coltemp<-c("#FF0000",OcolorZ[1:9])
coltemp[setdiff(2:(Nn+1),EdgN)]="#FFFFFF"  # Assign white to the nodes that are not connected
V(gN)$color=coltemp
coltemp[coltemp !="#FFFFFF"]="black"  # Assign white to the nodes that are not connected
V(gN)$frame.color=coltemp
nodenametemp<-as.character(1:(Nn+1))
nodenametemp[setdiff(2:(Nn+1),EdgN)]=""
V(gN)$name=nodenametemp
V(gN)$label.cex=3


# Signal only 
a=rep(1,length(EdgS))
gS <- graph(edges=c(rbind(a,EdgS)),n=10,directed=F)
gS <- gset(gS)
coltemp<-c("#FF0000",OcolorZ[10:18])
coltemp[setdiff(2:(Nn+1),EdgS)]="#FFFFFF"  # Assign white to the nodes that are not connected
V(gS)$color=coltemp
coltemp[coltemp !="#FFFFFF"]="black" 
V(gS)$frame.color=coltemp
nodenametemp<-as.character(1:(Nn+1))
nodenametemp[setdiff(2:(Nn+1),EdgS)]=""
V(gS)$name=nodenametemp
V(gS)$label.cex=3

# Plotting the pruned graph (After filtering case)
# Nominal
a=rep(1,length(EdgN_AF))
gNF <- graph(edges=c(rbind(a,EdgN_AF)),n=10,directed=F)
gNF <- gset(gNF)
coltemp<-c("#FF0000",OcolorZ[19:27])
coltemp[setdiff(2:(Nn+1),EdgN_AF)]="#FFFFFF"  # Assign white to the nodes that are not connected
V(gNF)$color=coltemp
coltemp[coltemp !="#FFFFFF"]="black"  # Assign white to the nodes that are not connected
V(gNF)$frame.color=coltemp
nodenametemp<-as.character(1:(Nn+1))
nodenametemp[setdiff(2:(Nn+1),EdgN_AF)]=""
V(gNF)$name=nodenametemp
V(gNF)$label.cex=3

# Filter only 
a=rep(1,length(EdgF))
gF <- graph(edges=c(rbind(a,EdgF)),n=10,directed=F)
gF <- gset(gF)
coltemp<-c("#FF0000",OcolorZ[28:36])
coltemp[setdiff(2:(Nn+1),EdgF)]="#FFFFFF"  # Assign white to the nodes that are not connected
V(gF)$color=coltemp
coltemp[coltemp !="#FFFFFF"]="black" 
V(gF)$frame.color=coltemp
nodenametemp<-as.character(1:(Nn+1))
nodenametemp[setdiff(2:(Nn+1),EdgF)]=""
V(gF)$name=nodenametemp
V(gF)$label.cex=3

# Filter + Signal
a=rep(1,length(EdgFS))
gFS <- graph(edges=c(rbind(a,EdgFS)),n=10,directed=F)
gFS <- gset(gFS)
coltemp<-c("#FF0000",OcolorZ[37:45])
coltemp[setdiff(2:(Nn+1),EdgFS)]="#FFFFFF"  # Assign white to the nodes that are not connected
V(gFS)$color=coltemp
coltemp[coltemp !="#FFFFFF"]="black" 
V(gFS)$frame.color=coltemp
nodenametemp<-as.character(1:(Nn+1))
nodenametemp[setdiff(2:(Nn+1),EdgFS)]=""
V(gFS)$name=nodenametemp
V(gFS)$label.cex=3
plot(gFS,layout=coord,mark.groups = list(c(3,4),c(8,9)),
     mark.col ="lightblue",vertex.label.family="serif",mark.border = F)

###### Fig 4b ##########

par(mar=c(0,0,0,0),mfrow=c(1,1))
plot(gN,layout=coord,mark.groups = list(c(3),c(5,7),c(8,9,10)),
     mark.col ="lightblue",vertex.label.family="serif",mark.border = F)
text(0,-.93,"Seed",cex = 3,col="blue",font=2)

###### Fig 4c ##########

par(mar=c(0,0,0,0),mfrow=c(1,1))
plot(gS,layout=coord,mark.groups = list(c(3),c(8,9)),
     mark.col ="lightblue",vertex.label.family="serif",mark.border = F)
text(0,-.93,"Seed",cex = 3,col="blue",font=2)

###### Fig 4e ##########

par(mar=c(0,0,0,0),mfrow=c(1,1))
plot(gNF,layout=coord,mark.groups = list(c(3,4),c(5,7),c(8,9,10)),
     mark.col ="lightblue",vertex.label.family="serif",mark.border = F)
text(0,-.93,"Seed",cex = 3,col="blue",font=2)

###### Fig 4f ##########

par(mar=c(0,0,0,0),mfrow=c(1,1))
plot(gF,layout=coord,mark.groups = list(c(3,4),c(5,7),c(8,10)),
     mark.col ="lightblue",vertex.label.family="serif",mark.border = F)
text(0,-.93,"Seed",cex = 3,col="blue",font=2)

###### Fig 4g ##########

par(mar=c(0,0,0,0),mfrow=c(1,1))
plot(gFS,layout=coord,mark.groups = list(c(3,4),c(8,9)),
     mark.col ="lightblue",vertex.label.family="serif",mark.border = F)
text(0,-.93,"Seed",cex = 3,col="blue",font=2)

###### Color bar for  Fig 4b,c,e,f,g ##########
ma=max(Zarr)
mi=min(Zarr)
scale = (length(OcolorZZ))/(ma-mi)
ticks=round(seq(mi,ma, len=11))
plot(20,20, type='n', xaxt='n', bty='n',xlab='Z score',cex.lab=2.2, yaxt='n', ylab='',xlim=c(mi,ma),ylim=c(0,60))
axis(side=1,at=ticks,cex.axis=2)

for (i in 1:(length(OcolorZZ))) {
  y = (i-1)/scale + mi
  rect(y,-3,y+1,2.5,col=OcolorZZ[i], border=NA)    # Change the last value here for the height of the color bar
}	
