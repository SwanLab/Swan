function [VAR,celasST]= SmallStrainJ2Plasticity(VAR,MATPRO,DATA,VARint_n)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
if nargin == 0
    load('tmp.mat')
end


% The function below is a copy UpdateStresses_J2_FE2D.m 
[celasST,  VAR.PK2STRESS, VAR.PlasticStrains,  VAR.InternalVarStrain,VAR.YieldStress ] = ...
    SmallStrainJ2Plasticity_LOCAL(VAR.GLSTRAINS,MATPRO,VARint_n.PlasticStrains,...
    VARint_n.YieldStress,VARint_n.InternalVarStrain,...
    DATA)  ; 

 