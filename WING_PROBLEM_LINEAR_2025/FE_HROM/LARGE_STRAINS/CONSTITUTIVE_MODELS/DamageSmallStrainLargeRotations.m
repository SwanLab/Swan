function [VAR,celasST]= DamageSmallStrainLargeRotations(VAR,STRAIN,MATPRO,VARINT)
%==========================================================================
% function [VAR, celasST] = DamageSmallStrainLargeRotations(VAR, STRAIN, MATPRO, VARINT)
%--------------------------------------------------------------------------
% PURPOSE:
%   Implements an **isotropic scalar damage model** under small strain but
%   allowing for **large rotations** at the global (finite element) level.
%   The model follows a standard energy-equivalence assumption:
%
%       σ = (1 - d) * σ_eff
%
%   where:
%       σ_eff  = effective (undamaged) stress tensor
%       d      = scalar damage variable,  0 ≤ d < 1
%
%   The algorithm performs an **incremental update** of the internal
%   variables:
%       - r : strain-like variable (damage driving force)
%       - q : stress-like variable (hardening variable)
%       - d : damage variable
%
%   and computes both the **current stress state** and the **consistent
%   tangent operator** (algorithmic stiffness) for integration at Gauss
%   points.
%
%--------------------------------------------------------------------------
% INPUTS:
%   VAR      : Structure containing the field VAR.PK2STRESS (initialized).
%   STRAIN   : Strain tensor at Gauss points (Voigt vector, possibly
%              vectorized for multiple points).
%   MATPRO   : Structure containing material properties:
%                 - MATPRO.celasglo : Elasticity tensor of the undamaged
%                                     material, in Voigt notation.
%                 - MATPRO.Hmodul   : Hardening modulus (per Gauss point).
%
%   VARINT   : Structure of internal variables at the previous time step:
%                 - VARINT.r_IntVARstrain : strain-like variable (r_n)
%                 - VARINT.q_IntVARstress : stress-like variable (q_n)
%                 - VARINT.d_DAMAGE       : damage variable (d_n)
%
%--------------------------------------------------------------------------
% OUTPUTS:
%   VAR      : Updated structure containing:
%                 - VAR.PK2STRESS : Second Piola–Kirchhoff stress tensor
%                                   at each Gauss point, including damage.
%
%   celasST  : Updated **algorithmic tangent matrix**, incorporating
%              damage effects and the additional term:
%
%                 C^load = (1 - d) * C
%                           + ((H * r - q) / r^3) * (σ_eff ⊗ σ_eff)
%
%              (see "DamageModels2025_main.pdf" for derivation)
%
%--------------------------------------------------------------------------
% ALGORITHM OUTLINE:
%   1. Compute the undamaged (effective) stress: σ_eff = C * ε.
%   2. Evaluate the damage driving variable r_trial.
%   3. Detect points in the elastic regime (r_trial ≤ r_n).
%      → Assign undamaged response, σ = (1 - d_n) * σ_eff.
%   4. For points entering the damaged regime:
%        - Update internal variables (r, q, d).
%        - Compute damaged stress: σ = (1 - d) * σ_eff.
%        - Build consistent tangent operator including the σ_eff⊗σ_eff term.
%
%--------------------------------------------------------------------------
% REFERENCES:
%   - J.A. Hernández Ortega, "Damage Models 2025" (internal notes),
%     /DEVELOPMENTS/DamageModels2025_main.pdf
%   - Simo, J.C. & Ju, J.W., "Strain- and stress-based continuum damage
%     models. I. Formulation", Int. J. Solids Struct. 23(7), 1987.
%
%--------------------------------------------------------------------------
% REMARKS:
%   • This formulation assumes small strains but large rotations at the
%     element level (the stress tensor is referred to the material frame).
%   • Fully vectorized implementation: all Gauss points are processed at once.
%   • The stress-tensor outer product (σ_eff⊗σ_eff) is computed via
%     the auxiliary function `DamageMatStressStress`.
%
%--------------------------------------------------------------------------
% HISTORY:
%   - 21-Oct-2025, J.A. Hernández Ortega
%     First implementation. Adapted from the isotropic linear hardening
%     model in /112_NonLIN_ROM_RBF/17_DAMAGE.mlx
%--------------------------------------------------------------------------
%==========================================================================


% Damage model, isotropic, linear hardening
% JAHO, 21-Oct-2025,
% SEe /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/17_DAMAGE.mlx
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/DamageModels2025_main.pdf
if nargin == 0
    load('tmp1.mat')
end

r_n = VARINT.r_IntVARstrain; % Strain-like internal variable
q_n = VARINT.q_IntVARstress; % Stress-like internal variable
d_n = VARINT.d_DAMAGE;       % damage

r_np1 = r_n;
q_np1 = q_n;
d_np1 = d_n;




STRESSeff = StressFromStrain_LINEAR(STRAIN,MATPRO) ; % Effective stress, undamaged material

nstrain = size(MATPRO.celasglo,2) ;

r_trial = SymmetricDamageModel(STRAIN,STRESSeff,nstrain) ;


% Points in the elastic regime
IndexPointsElastic = find(r_trial <= r_n) ;
IndexCompElastic = small2large(IndexPointsElastic,nstrain) ;
%VAR.PK2STRESS = zeros(size(STRESSeff)) ;
celasST  = MATPRO.celasglo ;
for istrain = 1:nstrain
    icomp= IndexCompElastic(istrain:nstrain:length(IndexCompElastic));
    VAR.PK2STRESS(icomp) = (1-d_n(IndexPointsElastic)).*STRESSeff(icomp) ;
    celasST(icomp,:) = bsxfun(@times,MATPRO.celasglo(icomp,:),(1-d_n(IndexPointsElastic)) ) ; 
end
%  VAR.r_IntVARstrain =r_n  ;
%  VAR.q_IntVARstress =q_n  ;
%  VAR.d_DAMAGE =d_n  ;


IndDAM = setdiff(1:length(r_trial),IndexPointsElastic) ;
if ~isempty(IndDAM)
    
    r_np1(IndDAM) = r_trial(IndDAM) ;
    q_np1(IndDAM) = q_n(IndDAM) + MATPRO.Hmodul(IndDAM).*(r_np1(IndDAM) - r_n(IndDAM)) ;
    d_np1(IndDAM) = 1-q_np1(IndDAM)./r_np1(IndDAM) ;
    IndDAM_comp = small2large(IndDAM,nstrain) ;
    
    
    
    
    % Stresses
    for istrain = 1:nstrain
        icomp= IndDAM_comp(istrain:nstrain:length(IndDAM_comp));
        VAR.PK2STRESS(icomp) = (1-d_np1(IndDAM)).*STRESSeff(icomp) ;
         celasST(icomp,:) = bsxfun(@times,MATPRO.celasglo(icomp,:),(1-d_np1(IndDAM)) ) ; 
    end    
    
    % ALGORITHMIC TANGENT MATRIX, component dependent on
    % sigmaEFF*sigmaEFF^T
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/DamageModels2025_main.pdf
    % \C^{load} = (1-\dDMG)  \C   +  \Par{\dfrac{\Hhs \rDMG - \qDMG}{\rDMG^3}}  \stressEFF   \stressEFF^T
    FactorSS = MATPRO.Hmodul(IndDAM).*r_np1(IndDAM) - q_np1(IndDAM);
    FactorSS = FactorSS./r_np1(IndDAM).^3; 
    MatStressStress = DamageMatStressStress(STRESSeff(IndDAM_comp),nstrain) ; 
     for istrain = 1:nstrain
         icompLOC = istrain:nstrain:length(IndDAM_comp) ; 
        icomp= IndDAM_comp(icompLOC);  
        
         celasST(icomp,:) =    celasST(icomp,:) + bsxfun(@times,MatStressStress(icompLOC,:),FactorSS) ; 
    end   
    
end

  VAR.r_IntVARstrain = r_np1; % Strain-like internal variable
  VAR.q_IntVARstress = q_np1; % Stress-like internal variable
  VAR.d_DAMAGE = d_np1;       % damage


