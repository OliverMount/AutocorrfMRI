# AutocorrfMRI
R codes to reproduce all the figures for the manuscript The Impact of sampling rate in variance correction for single subject fMRI connectivity analysis.



The programs are in the following order

Fig. 2 and 4. The program is used to produce the Fig 2 and subfigures in Fig. 4  [Refer Fig2_Fig4.R file]. 

Fig. 3a. Empirical PDF of the Z score for filtered and unfiltered data [Refer Fig3a.R]

Fig. 3b. Empirical PDF of the t score for filtered and unfiltered data [Refer Fig4a.R]

Fig. 5a. Variance of SPCC due to filter  for various TRs [AFNI (IIR) filters]  
(Please extract the Filter.zip and place them in the folder in the location of the Fig5a.R ) [Refer Fig5a.R]

Fig. 5b. Variance of SPCC due to filter for various TRs [AFNI (IIR), FSL and FIR] 
(Please extract the Filter.zip and place them in the folder in the location of the Fig5b.R ) [Refer Fig5b.R]


Fig. 6. Histograms of variance of the SPC before and after filtering 
( Please extract the NKI_DOF_BF.zip  and NKI_DOF_AF.zip and place the extracted folders ( NKI_DOF_BF and NKI_DOF_AF) in the  location same as  Fig6.R  or change the path accordingly in Fig6.R)  [Refer Fig6.R]

Fig. 7. False positive rate of variance correction methods [Refer Fig7.R]



Fig. 8. Seed-voxel-based connectivity in the DMN    [Refer Fig8.R]

Fig. 9. Seed-voxel-based connectivity in the DMN  (prospective and retrospective measurements)  [Refer Fig9.R]



VoxelDOF.R  -- to compute the variance (1/DOF) between a seed voxel and another voxel time series

   
 


The Files here will be updated soon
