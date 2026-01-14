function [PROPMAT,DATA ]=  MaterialDATA_elas3D_2mat(DATA)
% -----------------------------------------------------------------------------------
% 3. Material data  (linear elasticity)
% -----------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------
% 3. Material data  (linear elasticity)
% -----------------------------------------------------------------------------------
imat =1 ; % Index material
% Elasticity matrix
E = 70000  ; %  MPa, Young's modulus
nu = 0.3; % Poisson's coefficient
% Compliance matrix for an isotropic materials (with all entries, 3D)
% See slides, page 23.% figure(1)
% hold on
% xlabel('Time')
% ylabel('Amp disp')
% for iloadstate = 1:length(DIRICHLET(icond).PRESCRIBED_DISP)
%     fplot(DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).TIMEFUN,DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).INTERVAL)
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = E/2/(1+nu) ;  % Shear modulus
celasINV3D = [1/E  -nu/E -nu/E  0 0 0
    -nu/E  1/E  -nu/E  0 0 0
    -nu/E  -nu/E  1/E  0 0 0
    0      0    0  1/G   0 0
    0      0      0  0 1/G 0
    0      0      0  0  0  1/G] ;
ElasticityMatrix = inv(celasINV3D) ;
PROPMAT(imat).ElasticityMatrix =  ElasticityMatrix  ; %
PROPMAT(imat).Density = 2.7e-3 ; % MKg/m3
DATA.TYPE_CONSTITUTIVE_MODEL_ALL = 'SMALL_STRAINS_ELASTIC' ;  
%DATA.StrainStressWith4Components = 1; 

  
 IREF = 1;   % REference material (top plate) (hard, 700000 MPa Young's Modulus)
 
imat  =2 ;  % Core material, soft
%FACTOR_red = 1; 
PROPMAT(imat) = PROPMAT(IREF) ;  
 % PROPMAT(imat).ElasticityMatrix = FACTOR_red*ElasticityMatrix  ;  
  
% imat  =3 ;   % Similar to material 1 (hard)
% PROPMAT(imat) = PROPMAT(IREF) ;  
%   FACTOR_red = 1; 
%   PROPMAT(imat).ElasticityMatrix= FACTOR_red*ElasticityMatrix  ;  
  
% imat  =4 ;  % Top/bottom plates, hard
% PROPMAT(imat) = PROPMAT(IREF) ;  
%   FACTOR_red = 1; 
%   PROPMAT(imat).ElasticityMatrix = FACTOR_red*ElasticityMatrix  ;   
%   
%   
%   imat  =5 ;  % SOft,
% PROPMAT(imat) = PROPMAT(IREF) ;  
%   FACTOR_red = 0.1; 
%   PROPMAT(imat).ElasticityMatrix = FACTOR_red*ElasticityMatrix  ;  
%  
%  