function [DATA,INICOND,MATPRO] = DamageMaterialPropIntPoints(DATA,MESH,PROPMAT)
%==========================================================================
% DamageMaterialPropIntPoints
%
% PURPOSE:
%   Preprocess material data for a small-strain isotropic damage model at
%   **integration-point level** and set consistent initial conditions for the
%   internal variables. This routine:
%     (i)   builds element elasticity matrices and pointwise material props
%           via SmallStrainDamagePROP,
%     (ii)  expands/assembles the element elasticities into a global
%           block-diagonal elasticity operator at Gauss-point resolution, and
%     (iii) defines the list of internal variables and initializes them.
%
% INPUTS:
%   DATA    - Problem-level structure. Required fields:
%               • DATA.MESH.ngaus   : # Gauss points per element
%               • DATA.MESH.nstrain : # strain components (3 in 2D, 6 in 3D)
%               • DATA.typePROBLEM  : 'pstrain' | 'pstress' | '3D'
%   MESH    - Mesh structure. Required fields:
%               • COOR          : nodal coordinates (nnode × ndim)
%               • MaterialType  : (nelem × 1) integer material id per element
%             Optional:
%               • nstrain       : if present, overrides DATA.MESH.nstrain
%   PROPMAT - Array of material property structs (one per material id).
%             Required fields (used indirectly by SmallStrainDamagePROP):
%               • ElasticityMatrix (6×6 in Voigt 3D)
%               • Density
%               • Strength
%               • HmodulusHARD
%
% OUTPUTS:
%   DATA    - Same as input, with:
%               • DATA.ListFieldInternalVariables : cell array with the names
%                 of the internal variables used by the constitutive update:
%                   {'q_IntVARstress','r_IntVARstrain','d_DAMAGE'}
%   INICOND - Structure of initial conditions at Gauss points:
%               • INICOND.d_DAMAGE       : zeros(ngaus_total,1)  (damage = 0)
%               • INICOND.r_IntVARstrain : rDMG_0  (energy-like threshold)
%               • INICOND.q_IntVARstress : rDMG_0  (stress-like conjugate)
%             where ngaus_total = nelem * DATA.MESH.ngaus.
%   MATPRO  - Material property structure at Gauss-point level:
%               • celasglo : global block-diagonal elasticity (expanded)
%               • dens, Strength, Hmodul, rDMG_0, ...
%
% METHOD / STEPS:
%   1) Ensure consistency of the number of strain components:
%        MESH.nstrain = DATA.MESH.nstrain
%   2) Call SmallStrainDamagePROP to:
%        – map per-material properties to elements and their Gauss points,
%        – compute rDMG_0 from a unit uniaxial tensile energy norm,
%        – return per-element elasticities MATPRO.celasglo(:,:,e).
%   3) Build the **global Gauss-point elasticity operator** by expanding the
%      element-level blocks for each element and its DATA.MESH.ngaus points:
%        Cglo = DefineElastMatGLO_nw(MATPRO.celasglo, DATA.MESH.ngaus)
%      and overwrite MATPRO.celasglo with the expanded/global version.
%   4) Register the internal variable names in DATA.ListFieldInternalVariables.
%   5) Set initial conditions:
%        • d_DAMAGE = 0 (undamaged material),
%        • r_IntVARstrain = rDMG_0,
%        • q_IntVARstress = rDMG_0.
%      (These choices initialize the model at the elastic threshold so that
%       the first increment can activate/compare against the same scale.)
%
% CONVENTIONS / SHAPES:
%   - celasglo (global): typically stored as a block-diagonal operator whose
%     blocks are (nstrain×nstrain), one block per Gauss point.
%   - Vectors Strength, Hmodul, rDMG_0, d_DAMAGE, q_IntVARstress,
%     r_IntVARstrain are sized (ngaus_total × 1).
%
% DEPENDENCIES:
%   SmallStrainDamagePROP.m
%   DefineElastMatGLO_nw.m   % expands/assembles element elasticities to GP level
%
% SEE ALSO:
%   SmallStrainDamage_LOCAL.m   % constitutive update (return mapping/damage)
%   SmallStrainJ2PlasticityPROP.m
%
% AUTHOR / DATE:
%   Joaquín A. Hernández (UPC/CIMNE), Barcelona, 21-Oct-2025
%==========================================================================

if nargin == 0
    load('tmp1.mat')
end


% See input file (example):
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/17_DAMAGE/DamageMaterialDATA.m
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/17_DAMAGE.mlx


MESH.nstrain=  DATA.MESH.nstrain  ;
[MATPRO] = SmallStrainDamagePROP(MESH,DATA.typePROBLEM,PROPMAT,DATA) ;
disp('Global elasticity matrix    ...')
Cglo = DefineElastMatGLO_nw(MATPRO.celasglo,DATA.MESH.ngaus) ; % ; ...
%  Cglo = ConvertBlockDiag(Cglo) ;
MATPRO.celasglo = Cglo ;
DATA.ListFieldInternalVariables = {'q_IntVARstress','r_IntVARstrain','d_DAMAGE'} ;
INICOND.d_DAMAGE = zeros(size(MATPRO.Hmodul)) ;  % damage = 0, initially
 INICOND.r_IntVARstrain = MATPRO.rDMG_0 ; 
  INICOND.q_IntVARstress = MATPRO.rDMG_0 ; 
