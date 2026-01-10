clc
clear all

syms a   Sx Sy Sz Syz Sxz Sxy

% AROUND Y-AXIS
R = [cos(a), 0  , -sin(a)
     0       1     0
     sin(a)  0     cos(a)] ; 
% AROUND X-AXIS
R = [1  0   0
    0  cos(a), -sin(a)
     0 sin(a)       cos(a)] ; 
 
 
 sigmaLOC = [Sx Sxy Sxz
          Sxy Sy Syz 
          Sxz Syz Sz] ; 
      
sigmaGLO = R*sigmaLOC*R.' ; 

% Row 1 

sigma_x =  factor(sigmaGLO(1,1),[Sx,Sy,Sz,Sxz,Syz,Sxy])  
sigma_y =   factor(sigmaGLO(2,2),[Sx,Sy,Sz,Sxz,Syz,Sxy])  
sigma_z =  factor(sigmaGLO(3,3),[Sx,Sy,Sz,Sxz,Syz,Sxy])  
tau_yz =  factor(sigmaGLO(2,3),[Sx,Sy,Sz,Sxz,Syz,Sxy])  
tau_xz =  factor(sigmaGLO(1,3),[Sx,Sy,Sz,Sxz,Syz,Sxy])  
tau_xy =  factor(sigmaGLO(1,2),[Sx,Sy,Sz,Sxz,Syz,Sxy])  




      
      
 