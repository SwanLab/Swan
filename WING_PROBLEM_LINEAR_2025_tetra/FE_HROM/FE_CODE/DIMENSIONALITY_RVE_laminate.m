clc
clear all

% Compute Generalized elasticity matrix (LAMINATE)
% -------------------------

RECALCULATE_Celas = 1 ;
NameFileMesh ='lam10lay_8000.msh' ; 'layer10_16rve.msh'; 'layer10_4rve_32000.msh';'layer10_16rve.msh' ;'lam10lay_128.msh' ;'lam10lay_32000.msh' ;
NameWS_unitcell ='DATAWS/Celas_iSO.msh.mat' ;  'DATAWS/Celas_ejem1.msh.mat' ; 'DATAWS/Celas_MYFIRSTMESH.msh.mat' ;
% Elasticity matrix for each  layer (and angle of rotation)
% The name of the .mat containg such information is to be provided
% -----------------------
%
nameWSstrain =  ['DATAWS/EJEMPLO_paper_laminate.mat'] ;

nplies = 8 ;
ANG_FIB = [0 45 0 90 0 60 0 -45]
for iply = 1:nplies
    MATERIAL.PLY(iply).NAMEWS = NameWS_unitcell;  % Store the Celas in this folder
    MATERIAL.PLY(iply).ANGLE = ANG_FIB(iply) ;  % Angle subtended by the fiber and the x-axis
end
DATA.RECALCULATE_STIFFNESS = 1;   % To avoid computing again the stiffness matrix (if =0) when computing
% the average stresses for different input
% strains, yet similar
% material/geometric properties
DATA.CALCULATE_averageSTRESS = 2; %

DATA.BOUNDARY_CONDlam = 'ZERO'; 'ZERO_minTB' ;    %PERIODIC';     % PERIODIC/ZERO
DATA.niterCONJG = 3000 ; % Number of iterations  conjugated gradient
% ---- END INPUTS ----------------------------------
STRAIN_GLO = {} ;
STRESS_GLO = {} ;
DATA.CALCULATE_STRENGTH_LAM = 0  ;

for i = 1:8
    strainINP = zeros(8,1) ;
    strainINP(i) = 1;
    [stressMACRO DATAOUT DATA]= LaminateSTRESScal(strainINP,NameFileMesh,MATERIAL,DATA) ;
    Celas(:,i) = stressMACRO ;
    DATA.RECALCULATE_STIFFNESS =0 ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
    STRAIN_GLO{i} = DATAOUT.strain ;
    STRESS_GLO{i} = DATAOUT.stress ;
    d_GLO{i} = DATAOUT.d ;
 end

%densCOMP = DATAOUT.densCOMP ;
load(DATA.nameWORKSPACE,'Bst') ;
save(nameWSstrain,'STRAIN_GLO','STRESS_GLO','DATAOUT','d_GLO','Bst') ;

