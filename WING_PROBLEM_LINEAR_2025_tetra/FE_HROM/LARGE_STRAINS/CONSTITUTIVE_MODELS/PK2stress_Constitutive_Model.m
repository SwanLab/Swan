function [StwoST,celasST,detFgrad ]= PK2stress_Constitutive_Model(EgreenlST,MATPRO,DATA,FgradST)

if nargin == 0
    load('tmp2.mat')
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
            [StwoST,celasST]= SmallStrainLargeRotations(EgreenlST,MATPRO) ;
        case 'NeoHookean'
            [StwoST,celasST]= NeoHookStressStrain(detFgrad,EgreenlST,MATPRO,DATA) ;
            
        case 'SMALL_STRAINS_J2_PLASTICITY'
            [StwoST,celasST]= SmallStrainJ2Plasticity(EgreenlST,MATPRO) ;
        otherwise
            error('OPtion not implemented yet')
    end    
end
