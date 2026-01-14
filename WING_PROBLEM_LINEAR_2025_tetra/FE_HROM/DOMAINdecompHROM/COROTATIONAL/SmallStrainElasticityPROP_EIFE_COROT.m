function [MATPRO,DATA] = SmallStrainElasticityPROP_EIFE_COROT(MESH,typePROBLEM,PROPMAT,DATA)
%--------------------------------------------------------------------------
%  SmallStrainElasticityPROP_EIFE_COROT
%
%  Computes the **elemental linear elasticity matrices** and **mass density vectors**
%  for each coarse-scale EIFEM element in the context of the **co-rotational formulation**
%  under **small strain** assumptions. This function is adapted from the standard version
%  to accommodate the **CECM (Continuous Empirical Cubature Method)**, which permits
%  a variable number of integration points per element.
%
%  INPUTS:
%    - MESH      : structure containing geometry and material type for each element
%    - typePROBLEM : string ('3D', 'pstrain', etc.)
%    - PROPMAT   : array of material data, each with EIFE_prop including CECM weights
%    - DATA      : user-defined input structure, including CECM integration point data
%
%  OUTPUTS:
%    - MATPRO : structure containing:
%         > celasglo : global block-diagonal elasticity matrix 
%                      (stacked per element and integration point)
%         > dens     : global density vector, aligned with RHS Gauss points
%    - DATA   : possibly enriched with defaults
%
%  THEORY:
%    - The elasticity matrices are assembled based on the number of **CECM internal force
%      integration points** defined for each EIF parent domain.
%
%    - The linear elasticity matrices are taken from the pre-trained library via:
%
%         celasglo_elem ← AssignOneEIFelementELASTICITYmatCOROT(...)
%
%    - This structure allows for **material heterogeneity across elements**, supporting
%      multi-type element simulations.
%
%  CONNECTION TO THEORY:
%    - The elasticity matrix corresponds to **Cf0′** in Eq. (404) of *EIFEM_largeROTfinal.pdf*:contentReference[oaicite:0]{index=0}
%
%    - The full stiffness matrix is constructed in:
%
%         KΦ,lin = λ^(d−2) * KΦ,lin_ref      (Eq. 409):contentReference[oaicite:1]{index=1}
%
%    - The integration weights and points are determined offline by the **CECM**, and
%      the actual B-matrix and elasticity matrix are applied during preprocessing.
%
%  CECM ROLE:
%    - This version assumes that the **number and location of integration points** can
%      vary element-wise, as opposed to the standard uniform Gaussian quadrature.
%    - It uses **per-element weighting** for both stiffness and mass matrix evaluation.
%
%  AUTHOR:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    Version: 28-Oct-2024, UPC Terrassa
%    Comments by ChatGPT4, 13-May-2025
%
%  SEE ALSO:
%    - AssignOneEIFelementELASTICITYmatCOROT
%    - AssignOneEIFelementDENSmat
%    - Section 8.3 and Eq. (409–410) of *EIFEM_largeROTfinal.pdf*:contentReference[oaicite:2]{index=2}
%
%--------------------------------------------------------------------------




% Adaptation of  SmallStrainElasticityPROP_EIFE for the Empirical Interscale FE
% method using the corotational approach (essentially the same as this function, but now able to handle variable
% number CECM integration points for internal forces)
% JAHO, 28-Oct-2024, Monday, UPC Terrassa.
% -----------------------------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end

% See also function DefineElastMatGLO_nw

ndim = size(MESH.COOR,2)  ; % Number of spatial dimensions

% This nstrain refers to the number of Cauchy stress components
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

% LINEAR ELASTICITY MATRIX
%

nelem = size(MESH.MaterialType,1) ;
%number_cecmpoints_INTERNAL_FORCES_all = sum(DATA.MESH.number_cecmpoints_INTERNAL_FORCES) ;
MATPRO.celasglo = cell(nelem,1) ;    % zeros(nstrain*number_cecmpoints_INTERNAL_FORCES_all,nstrain) ;  % Global array of elasticity matrices
MATPRO.dens = zeros(nelem*DATA.MESH.ngaus_RHS,1) ;
%celasgloINV = zeros(6,6,nelem) ;
for imat = 1:length(PROPMAT)
    % Elasticity matrix for one type of EIF element
    % ---------------------------------------------------------------------------------------------------------
    PROPMATLOC = PROPMAT(imat).PROPMAT ;
    MaterialTypeLocal = PROPMAT(imat).EIFE_prop.INTforces.MaterialType ;
   % ncecmINT_forces =  length(PROPMAT(imat).EIFE_prop.INTforces.weights) ;
    
    celasglo_elem = AssignOneEIFelementELASTICITYmatCOROT(nstrain,DATA,PROPMATLOC,MaterialTypeLocal,typePROBLEM)  ;
    %----------------------------------------------------------------------------------------------------------
    
    ELEMS = find(MESH.MaterialType == imat) ; % Elements of this type of material
    
    MATPRO.celasglo(ELEMS) = {celasglo_elem} ;
    
    
    if ~isempty(MESH.VOLUME)
    MaterialTypeLocal = PROPMAT(imat).EIFE_prop.BodyForces.MaterialType ;
    dens_elem = AssignOneEIFelementDENSmat(DATA,PROPMATLOC,MaterialTypeLocal)  ;
    ELEMSdofsDENS = small2large(ELEMS,size(dens_elem,1)) ;
    MATPRO.dens(ELEMSdofsDENS,:) = repmat(dens_elem,length(ELEMS),1)  ;
    end
    
    
end

MATPRO.celasglo = cell2mat(MATPRO.celasglo) ;

