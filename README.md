# AutocorrfMRI

Autocorrelation both inherent and pre-processing induced in the fMRI time series affects the connectivity analysis. The autocorrelation issue is largely ignored in the fMRI preprocessing literature. We bring to forefront the impact of the autocorrelation in connectivity analysis using real and simulated data.


This repository contains R programs to reproduce the figures in the paper: "The Impact of sampling rate on statistical significance for single subject fMRI connectivity analysis."

https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.24600

----------------------------------------------------------------------------------------

Fig2andFig4.R  ----    To produce Fig 2 and subfigures in Fig. 4
 

Fig3a.R   ----    To reproduce empirical PDF a) Z-score 

Fig3b.R   ----    To reproduce empirical PDF b) t-score

Fig5a.R   ----     To produce variance of SPCC due to filter for various TRs [AFNI (IIR) filters]  
                (Please extract the Filter.zip and place them in the folder in the location of the Fig5a.R ) 

Fig5b.R   ----     To produce variance of SPCC for AFNI (IIR), FSL and FIR filters  
                (Please extract the Filter.zip and place them in the folder in the location of the Fig5b.R ) 


Fig6.R ----  Histograms of variance of the SPC before and after filtering 
(Please extract the NKI_DOF_BF.zip  and NKI_DOF_AF.zip and place the extracted folders ( NKI_DOF_BF and NKI_DOF_AF) in the  location same as  Fig6.R  or change the path accordingly in Fig6.R)   

Fig7.R ---- False positive rate of variance correction methods [Refer Fig7.R]



Fig8.R  ---- Seed-voxel-based connectivity in the DMN    [Refer Fig8.R]

Fig9.R  ----  Seed-voxel-based connectivity in the DMN  (prospective and retrospective measurements)  [Refer Fig9.R]

----------------------------------------------------------------------------------------
Additional function files
----------------------------------------------------------------------------------------

VoxelDOF.R  ----  to compute the variance (1/DOF) between a seed voxel and another voxel time series

VarSPCC.R  --- function to compute the variance of the SPCC between two time series x1 and x2.

(The files are periodically updated)
