function [MATPRO] = SmallStrainDamagePROP(MESH,typePROBLEM,PROPMAT,DATA)
%==========================================================================
% SmallStrainDamagePROP
%
% PURPOSE:
%   Builds the structure `MATPRO` containing all material properties 
%   and precomputed quantities required for a small-strain isotropic 
%   **damage model**. This function processes the material data defined 
%   in `PROPMAT` and assigns them to each finite element and Gauss point.
%
%   The routine also computes the reference energy norm associated with 
%   a unit uniaxial tensile stress state, used to define the initial 
%   damage threshold `rDMG_0`. This normalization ensures that the scalar 
%   internal variable `r` (strain-like) has consistent energy units.
%
% INPUTS:
%   MESH         - Structure containing mesh data:
%                  • COOR : nodal coordinates (nnode × ndim)
%                  • MaterialType : (nelem × 1) array mapping each element 
%                    to its corresponding material index in PROPMAT.
%                  • nstrain (optional) : number of strain components 
%                    (3 for 2D, 6 for 3D problems)
%
%   typePROBLEM  - String defining the mechanical formulation:
%                  'pstrain' : plane strain
%                  'pstress' : plane stress
%                  '3D'      : full 3D solid
%
%   PROPMAT      - Array of material property structures, each containing:
%                  • ElasticityMatrix : Elastic constitutive tensor (6×6)
%                  • Density           : Material density
%                  • Strength          : Tensile strength (stress-like)
%                  • HmodulusHARD      : Hardening modulus (0 < H < 1)
%
%   DATA         - Structure with additional problem-level data, 
%                  including number of Gauss points per element:
%                  • MESH.ngaus : number of Gauss points per element
%
% OUTPUT:
%   MATPRO       - Structure containing all processed material data:
%                  • celasglo : (nstrain × nstrain × nelem) array of 
%                    elasticity matrices for each element
%                  • dens     : (nelem × 1) density array
%                  • Strength : (ngaus × 1) strength at each Gauss point
%                  • Hmodul   : (ngaus × 1) hardening modulus
%                  • rDMG_0   : (ngaus × 1) reference damage threshold, 
%                    computed as Strength × EnergyNorm
%
% ALGORITHM NOTES:
%   1. For each material, the 3D elasticity matrix is reduced according 
%      to the selected mechanical assumption ('pstrain', 'pstress', '3D').
%   2. The energy norm is computed as:
%         EnergyNorm = sqrt(σᵗ·ε)
%      for a uniaxial tensile stress σ = [1, 0, 0, 0, 0, 0]ᵗ.
%   3. The resulting `rDMG_0` defines the scale relating stress-like 
%      and strain-like internal variables in the damage model.
%
% SEE ALSO:
%   SmallStrainDamage_LOCAL.m
%   SmallStrainJ2PlasticityPROP.m
%
% AUTHOR:
%   Joaquín A. Hernández, UPC/CIMNE
%   Barcelona, 21-Oct-2025
%==========================================================================

if nargin == 0
    load('tmp1.mat')
end

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
MATPRO.celasglo = zeros(nstrain,nstrain,nelem) ;  % Global array of elasticity matrices
MATPRO.dens = zeros(nelem,1) ;



ngausE = DATA.MESH.ngaus ;
ngaus =  nelem*ngausE ;
% E_young = zeros(ngaus,1) ;  % Young's modulus
% Hmodulus = zeros(ngaus,1) ; % Hardening parameters
% PoissonCOEF = zeros(ngaus,1) ;  % Poisson's ratio
% sigmay_0 = 1e20*ones(ngaus,1) ;  % Yield stress

MATPRO.Strength = zeros(ngaus,1) ;
MATPRO.Hmodul = zeros(ngaus,1) ;
MATPRO.rDMG_0 = zeros(ngaus,1) ;


%celasgloINV = zeros(6,6,nelem) ;
for imat = 1:length(PROPMAT)
    
    TENSILE_STRESS_3D  = zeros(6,1) ; 
TENSILE_STRESS_3D(1) = PROPMAT(imat).Strength ; 
    
    celas3D =PROPMAT(imat).ElasticityMatrix ; %
    INVcelas3D = inv(celas3D) ;
    TENSILE_STRAIN_3D = INVcelas3D*TENSILE_STRESS_3D ; 
    EnergyNorm = sqrt(TENSILE_STRESS_3D'*TENSILE_STRAIN_3D) ; 
    ELEMS = find(MESH.MaterialType == imat) ;
    switch typePROBLEM
        case 'pstrain'
            if nstrain == 3
                rowcol = [1 2 6] ;
            elseif nstrain == 4
                rowcol = [1 2 6 3] ;
            else
                error('Option not implemented')
            end
            celas = celas3D(rowcol,rowcol) ;
            
        case 'pstress'
            if nstrain == 3
                rowcol = [1 2 6] ;
                celasINV3D = inv(celas3D) ;
                celasINV = celasINV3D(rowcol,rowcol) ;
                celas = inv(celasINV) ;
            else
                error('Option not implemented')
            end
        case '3D'
            celas = celas3D ;
    end
    
    
    
    for eLOC=1:length(ELEMS)
        e = ELEMS(eLOC) ;
        pointsMAT = small2large(e,ngausE) ;
        
        MATPRO.celasglo(:,:,e) = celas ;
        MATPRO.dens(e) = PROPMAT(imat).Density ;
        MATPRO.Hmodul(pointsMAT) = PROPMAT(imat).HmodulusHARD  ;
        MATPRO.Strength(pointsMAT) = PROPMAT(imat).Strength  ;
         MATPRO.rDMG_0(pointsMAT) =  EnergyNorm  ;
    end
end
