IJK2MNI <- function(a,b)
{
  # Function to compute IJK (Voxel) coordinates to MNI (mm) 
  
  A=rbind(matrix(a,c(3,4),byrow=T),c(0,0,0,1))
  d=c(b,1)  
  VoxCor=round(A %*% d)  # A*d
  VoxCor[1:2,1]=-VoxCor[1:2,1]
  return(VoxCor[1:3]-abs(A[1,1]))    # AFNI order is different 
  
}