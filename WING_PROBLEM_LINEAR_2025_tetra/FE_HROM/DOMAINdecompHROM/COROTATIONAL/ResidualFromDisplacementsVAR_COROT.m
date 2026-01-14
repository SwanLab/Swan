function  [VAR,celastST,FgradST,detFgrad] =  ResidualFromDisplacementsVAR_COROT(OPERFE,VAR,MATPRO,DATA,VARint_n,...
    dCqLOC,KcLINq,BstQ,D_QrotALL) ;
% Adaptation of  ResidualFromDisplacementsVAR.m to the Co-rotational
% approach
% % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
% JAHO, 29-Oct-2024, Green's Aribau, Barcelona.  


if nargin == 0
    load('tmp3.mat')
end

[VAR,celastST,FgradST,detFgrad] = StressesFromDisplacementsVAR_COROT(OPERFE,VAR,MATPRO,DATA,VARint_n,...
    dCqLOC,BstQ,D_QrotALL) ;


if isempty(VAR.PK2STRESS)
    VAR.PoneST = [] ; Fint = [] ; VAR.RESID = [] ;
else    
    VAR.FINT = InternalForces_COROT(OPERFE,VAR.PK1STRESS,VAR.PK2STRESS,DATA,VAR,dCqLOC,KcLINq,BstQ,D_QrotALL) ;    
end
    % 6.2. Residual
    VAR.RESID  = VAR.FINT- VAR.FEXT;
end

