function  [VAR,celastST,FgradST,detFgrad] =  ResidualFromDisplacementsCABLES(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
%[GLSTRAINS,PK2STRESS,celastST,RESID,Fint,FgradST,PoneST,detFgrad] = ResidualFromDisplacements(OPERFE,d,MATPRO,DATA,FEXT)
% This is a copy of ResidualFromDisplacements.m.
% [VAR,celastST,Fint,FgradST,detFgrad] =...
%        ResidualFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA) ;

%The inputs are not
% specified explicitily, nor the outputs.
% See  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx


if nargin == 0
    load('tmp.mat')
end

% Stresses from displacements
%[GLSTRAINS,PK2STRESS,celastST,FgradST,PoneST,detFgrad] = StressesFromDisplacements(OPERFE,d,MATPRO,DATA) ;
%if DATA.INTERNAL_FORCES_USING_precomputed_Kstiff ==0
DATA.CALC_CTANG = 1; 
[VAR,celastST,FgradST,detFgrad] = StressesFromDisplacementsCABLES(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
 
VAR.FINT = InternalForcesCABLE(OPERFE,VAR.TENSIONV,DATA) ;

VAR.RESID  = VAR.FINT- VAR.FEXT;


  


