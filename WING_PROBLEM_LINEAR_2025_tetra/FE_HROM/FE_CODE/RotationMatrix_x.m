function T = RotationMatrix_x(R)





cos_a = R(2,2) ;
sin_a = R(3,2) ;

cosTc = (cos_a)^2 ;
sinTc = (sin_a)^2 ;
sinTd = 2*cos_a*sin_a ; %  (sin(2*a)) ;
cosTd = cosTc -sinTc ; % cos(2*a) ;


cosT = cos_a ;
sinT = sin_a ;

T = [cosTc  sinTc    0      0      0      -sinTd
    sinTc  cosTc    0      0      0      sinTd
    0      0       1      0      0        0
    0      0     0      cosT      sinT           0
    0      0        0       -sinT      cosT          0
    0.5*(sinTd)     -0.5*(sinTd)        0       0    0          cosTd         ] ;

% The above matrix is for a rotation around the z-axis. The desired
% rotation here is around the x-axis.   Hence we have to permute the
% indices so that 
%  OLD --> x,y,z
%  NEW --> y,z,x
% Therefore 
% OLD -->  xx, yy ,zz , yz, xz, xy
% NEW -->  zz, xx, yy , xy, yz, xz 
% ********  
PERM = [3,1,2,6,4,5] ; 
T = T(PERM,PERM) ; 

% This can be also inferred using symbolic function 
% /home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/FE_CODE/SymbolicStressRotation.m
