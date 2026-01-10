clc
clear all

% Compute elasticity matrix
% -------------------------

RECALCULATE_Celas = 1 ;
NameFileMesh = 'ejem1.msh'; 'ejemplo1.msh' ; % %
nameWSstrain =  ['DATAWS/EJEMPLO_paper2iso.mat'] ;
DATA.name_file_boundaryconditions = 'INPUT_zero'; 

% Transverse isotropic material
E = 230e3 ; nu = 0.3 ;
G = E/2/(1+nu)  ; 
MATERIAL.FIBER.E(1) = E ; %MPa  % Longitudinal elastic modulus
MATERIAL.FIBER.E(2) = E ; %MPa   % Transverse elastic modulus
MATERIAL.FIBER.nu(1) = nu ;  % Major poisson's ratio (contraction in the transverse directionupon an extension in the fiber direction)
MATERIAL.FIBER.nu(2) = nu ;  % Transverse poisson's ratio
MATERIAL.FIBER.Gshear(1) = G ; %MPa % In-plane shear modulus
MATERIAL.FIBER.Gshear(2)= G ; % MPa % Transverse shear modulus
MATERIAL.FIBER.INDEX = 2;  % Material index (within GID mesh)
MATERIAL.FIBER.strength  = 3930;  %  % Strength in MPa (same in all directions)
MATERIAL.FIBER.CONSTITUTIVE_MODEL  = 'VONMISES';  % s
MATERIAL.FIBER.DENSITY  =  1.8e-3 ; % MKg/m3 = 1e6 kg/m3 = 1e6 1e-3g/(1e6 cm3) = 1e-3 g/cm3    ;  % s
%%%
E = 3e3 ; 
MATERIAL.MATRIX.E(1) = E  ; % MPa Elastic modulus   % MATRIX
MATERIAL.MATRIX.nu(1) = nu ; % Poisson's ratio
MATERIAL.MATRIX.Gshear(1) =G ; % E/2/(1+nu)  ;  % In-plane shear modulus
MATERIAL.MATRIX.INDEX = 1; % Material index (within GID mesh)
MATERIAL.MATRIX.strength  = 64;  %   % Strength in MPa (same in all directions)
MATERIAL.MATRIX.CONSTITUTIVE_MODEL  = 'VONMISES';  %
MATERIAL.MATRIX.DENSITY  = 1.22e-3 ; % Mkg/m3  % s
%
% % Transverse isotropic material
% MATERIAL.FIBER.E(1) = 230e3 ; %MPa  % Longitudinal elastic modulus
% MATERIAL.FIBER.E(2) = 8e3 ; %MPa   % Transverse elastic modulus
% MATERIAL.FIBER.nu(1) = 0.25 ;  % Poisson's ratio (contraction in the transverse directionupon an extension in the fiber direction)
% MATERIAL.FIBER.nu(2) = 0.3 ;  % Transverse poisson's ratio
% MATERIAL.FIBER.Gshear(1) = 27.3e3 ; %MPa % In-plane shear modulus  G_LT
% MATERIAL.FIBER.Gshear(2)= 3.08e3 ; % MPa % Transverse shear modulus G_TT
% MATERIAL.FIBER.INDEX = 1;  % Material index
% %%%
% E = 3e3 ; nu = 0.3 ;
% MATERIAL.MATRIX.E(1) = E  ; % MPa Elastic modulus   % MATRIX
% MATERIAL.MATRIX.nu(1) = nu ; % Poisson's ratio
% MATERIAL.MATRIX.Gshear(1) = E/2/(1+nu)  ;  % In-plane shear modulus
% MATERIAL.MATRIX.INDEX = 2; % Material index


DATAdef.CALCULATE_AVGDENS = 1 ;
DATA.RECALCULATE_STIFFNESS =1 ;  % To avoid computing again the stiffness matrix (if =0) when computing
% the average stresses for different input
% strains, yet similar
% material/geometric properties
% ---- END INPUTS ----------------------------------
nameCELAS = ['DATAWS/Celas_',NameFileMesh,'.mat'] ;

%%%%%%%%%%%%%%%%
ntraj = 21 ;
TRAYECTORIES = zeros(6,ntraj) ;
TRAYECTORIES(:,1:6) = eye(6) ; 
Comb = nchoosek(1:6,2) ; 
for iii = 1:size(Comb,1)
    rowss = Comb(iii,:) ;
    TRAYECTORIES(rowss,iii+6) = 1 ; 
end
%  nad= 20 ; 
% TRAJ_AD = zeros(6,nad) ;
% Comb = nchoosek(1:6,3) ; 
% for iii = 1:size(Comb,1)
%     rowss = Comb(iii,:) ;
%     TRAJ_AD(rowss,iii) = 1 ; 
% end
% TRAYECTORIES = [TRAYECTORIES TRAJ_AD] ;
%%%%%%%%%%%%%%%
STRAIN_GLO = {} ;
STRESS_GLO = {} ;
d_GLO = {} ;

Celas = zeros(6,6) ;
for i = 1:size(TRAYECTORIES,2)
    %         strainINP = zeros(6,1) ;
    %         strainINP(i) = 1;
    strainINP = TRAYECTORIES(:,i) ;
    [ stressMACRO DATAOUT DATA]= CompHomog_CF(strainINP,NameFileMesh,MATERIAL,DATA) ;
    Celas(:,i) = stressMACRO ;
    DATA.RECALCULATE_STIFFNESS =0 ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
    STRAIN_GLO{i} = DATAOUT.strain ;
    STRESS_GLO{i} = DATAOUT.stress ;
    d_GLO{i} = DATAOUT.d ;
end
%densCOMP = DATAOUT.densCOMP ;
save(nameWSstrain,'STRAIN_GLO','STRESS_GLO','DATAOUT','d_GLO') ;
 