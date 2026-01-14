function [VAR,celasST,detFgrad ]= PK2stress_Constitutive_ModelVARincre(VAR,MATPRO,DATA,FgradST,VARint_n)
% Copy of  PK2stress_Constitutive_ModelVAR.m 
% % Adapted to determine "incremental" stresses (Nonlinear)
%  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
% 6-Jan-2024, Balmes, Barcelona 
%--------------------------------------

if nargin == 0
    load('tmp.mat')
elseif nargin == 3
    FgradST = [] ;
end
detFgrad = [] ;

% Determinant Gradient FgradST
if ~isempty(FgradST)
    detFgrad = Determinant_Fgrad(FgradST,DATA.MESH.ndim) ;
    
    IND_NEG = find(detFgrad<=0) ;
else
    IND_NEG = [] ;
end


%IND_NEG = [] ;
if ~isempty(IND_NEG)
    warning('Negative Jacobian...')
    IND_ELEM = large2smallREP(IND_NEG,DATA.MESH.ngaus) ;
    disp(['Elements with Negative Jacobians'])
    disp(IND_ELEM')
    %   clipboard('copy',num2str(IND_ELEM'))
    
    StwoST = [] ; celasST = [] ; detFgrad = [] ;
    
    
else
   
    switch  DATA.TYPE_CONSTITUTIVE_MODEL_ALL
        case 'SMALL_STRAINS_ELASTIC'
            [VAR.PK2STRESS,celasST]= SmallStrainLargeRotations(VAR.GLSTRAINS,MATPRO) ;
        case 'NeoHookean'
            [VAR.PK2STRESS,celasST]= NeoHookStressStrain(detFgrad,VAR.GLSTRAINS,MATPRO,DATA) ;
            
        case 'SMALL_STRAINS_J2_PLASTICITY'
            [VAR,celasST]= SmallStrainJ2Plasticity(VAR,MATPRO,DATA,VARint_n) ;
            
            
        otherwise
            error('OPtion not implemented yet')
    end    
end
