# AutocorrfMRI
This repository contains R programs to reproduce the figures in the manuscript: "The Impact of sampling rate in variance correction for single subject fMRI connectivity analysis."
----------------------------------------------------------------------------------------

Fig2andFig4.R  ----    To produce Fig 2 and subfigures in Fig. 4
 

Fig. 3a.R   ----    To reproduce empirical PDF a) Z-score 

Fig. 3b.R   ----    To reproduce empirical PDF b) t-score

Fig. 5a.R   ----     To produce variance of SPCC due to filter for various TRs [AFNI (IIR) filters]  
                (Please extract the Filter.zip and place them in the folder in the location of the Fig5a.R ) 

Fig. 5b.R   ----     To produce variance of SPCC for AFNI (IIR), FSL and FIR filters  
                (Please extract the Filter.zip and place them in the folder in the location of the Fig5b.R ) 


Fig. 6. ----  Histograms of variance of the SPC before and after filtering 
( Please extract the NKI_DOF_BF.zip  and NKI_DOF_AF.zip and place the extracted folders ( NKI_DOF_BF and NKI_DOF_AF) in the  location same as  Fig6.R  or change the path accordingly in Fig6.R)  [Refer Fig6.R]

Fig. 7. ---- False positive rate of variance correction methods [Refer Fig7.R]



Fig. 8. ---- Seed-voxel-based connectivity in the DMN    [Refer Fig8.R]

Fig. 9.----  Seed-voxel-based connectivity in the DMN  (prospective and retrospective measurements)  [Refer Fig9.R]



VoxelDOF.R  ----  to compute the variance (1/DOF) between a seed voxel and another voxel time series


(The files are periodically updated)
