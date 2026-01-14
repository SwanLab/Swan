function [T ]= RotationMatrixStress_y(R)

cos_a = R(1,1) ;
sin_a = R(3,1) ;
 
cosTc = cos_a^2 ; 
sinTc = sin_a^2 ; 
sin_2a = 2*cos_a*sin_a ; %  (sin(2*a)) ; 
cos_2a = cosTc -sinTc ; % cos(2*a) ; 



% See /home/joaquin/Desktop/CURRENT_TASKS/COMPOSITE_MATERIALS_DOCENCIA/APUNTES_LATEX/...
% DOCUMENTOS_anexos/MATLAB/ELASTOSTATIC_GEN/FE_CODE/SymbolicStressRotation.m

T  =zeros(6,6) ; 
%  Sx, Sy, Sz, Syz, Sxz, Sxy
%  1   2   3    4    5    6

% 1) Sx*cos(a)^2 + Sz*sin(a)^2 - 2*Sxz*cos(a)*sin(a)
i = 1; 
T(i,1) = cos_a^2  ;  T(i,3) = sin_a^2 ;  T(i,5) = -2*cos_a*sin_a;
% 2) Sy
i = 2; 
T(i,2) = 1  ;   
% 3) Sz*cos(a)^2 + Sx*sin(a)^2 + 2*Sxz*cos(a)*sin(a)
i = 3; 
T(i,3) = cos_a^2  ;  T(i,1) = sin_a^2 ;  T(i,5) = 2*cos_a*sin_a ; 
% 4) Syz*cos(a) + Sxy*sin(a)
i = 4; 
T(i,4) = cos_a  ;  T(i,6) = sin_a ;   
% 5)Sxz*cos(2*a) + (Sx*sin(2*a))/2 - (Sz*sin(2*a))/2
i = 5; 
T(i,5) = cos_2a  ;    T(i,1) = sin_2a/2  ;  T(i,3) = -sin_2a/2  ;
% 6)  Sxy*cos(a) - Syz*sin(a)
i = 6; 
T(i,6) = cos_a  ;  T(i,4) = -sin_a ;  %  





% 
% Tz = [cosTc  sinTc    0      0      0      -sinTd  
%   sinTc  cosTc    0      0      0      sinTd   
%      0      0       1      0      0        0       
%      0      0     0      cosT      sinT           0        
%      0      0        0       -sinT      cosT          0         
%     0.5*(sinTd)     -0.5*(sinTd)        0       0    0          cosTd         ] ;
% 
% % The above is a rotation about z. The expression for the rotationa around
% % y reads 
% Ty = zeros(size(Tz)) ;
