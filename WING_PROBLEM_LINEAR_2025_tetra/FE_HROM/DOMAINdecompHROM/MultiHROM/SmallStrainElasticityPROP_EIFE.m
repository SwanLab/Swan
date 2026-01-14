function [MATPRO,DATA] = SmallStrainElasticityPROP_EIFE(MESH,typePROBLEM,PROPMAT,DATA)
% Adaptation of  SmallStrainElasticityPROP for the Empirical Interscale FE
% method.
% JAHO, 21-March-2023
% -----------------------------------------
if nargin == 0
    load('tmp.mat')
end

% See also function DefineElastMatGLO_nw

ndim = size(MESH.COOR,2)  ;

MESH = DefaultField(MESH,'nstrain',[]) ;
if isempty(MESH.nstrain)
    if ndim==2
        nstrain = 3;
    else
        nstrain = 6 ;
        typePROBLEM ='3D' ;
    end
else
    nstrain = MESH.nstrain ;
end

nelem = size(MESH.MaterialType,1) ;
MATPRO.celasglo = zeros(nstrain*nelem*DATA.MESH.ngaus_STRESS,nstrain) ;  % Global array of elasticity matrices
MATPRO.dens = zeros(nelem*DATA.MESH.ngaus_RHS,1) ;
%celasgloINV = zeros(6,6,nelem) ;
for imat = 1:length(PROPMAT)
    % Elasticity matrix for one type of EIF element
    % ---------------------------------------------------------------------------------------------------------
    PROPMATLOC = PROPMAT(imat).PROPMAT ;
    MaterialTypeLocal = PROPMAT(imat).EIFE_prop.INTforces.MaterialType ;
    
    celasglo_elem = AssignOneEIFelementELASTICITYmat(nstrain,DATA,PROPMATLOC,MaterialTypeLocal,typePROBLEM)  ;
    %----------------------------------------------------------------------------------------------------------
    
    ELEMS = find(MESH.MaterialType == imat) ; % Elements of this type of material
    ELEMSdofs = small2large(ELEMS,size(celasglo_elem,1)) ;
    MATPRO.celasglo(ELEMSdofs,:) = repmat(celasglo_elem,length(ELEMS),1)  ;
    
    
    MaterialTypeLocal = PROPMAT(imat).EIFE_prop.BodyForces.MaterialType ;    
    dens_elem = AssignOneEIFelementDENSmat(DATA,PROPMATLOC,MaterialTypeLocal)  ;  
    ELEMSdofsDENS = small2large(ELEMS,size(dens_elem,1)) ;
    MATPRO.dens(ELEMSdofsDENS,:) = repmat(dens_elem,length(ELEMS),1)  ;
    
    
end



% Cglo = DefineElastMatGLO_nw(MATPRO.celasglo,DATA.MESH.ngaus) ; % ; ...
