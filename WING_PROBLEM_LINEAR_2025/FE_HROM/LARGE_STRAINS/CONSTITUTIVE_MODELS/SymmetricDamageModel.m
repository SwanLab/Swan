function r = SymmetricDamageModel(STRAIN,STRESSeff,nstrain)
%==========================================================================
% function r = SymmetricDamageModel(STRAIN, STRESSeff, nstrain)
%--------------------------------------------------------------------------
% PURPOSE:
%   Computes the **strain-like internal variable** `r` that drives the
%   evolution of isotropic damage in small-strain models.  
%
%   The variable `r` is defined as the square root of the scalar product
%   between the current strain tensor and the effective (undamaged) stress:
%
%       r = √( ε : σ_eff )
%
%   where:
%       ε       – strain tensor (Voigt notation)
%       σ_eff   – effective stress (undamaged material)
%
%   This symmetric energy-like definition ensures that `r` has units of
%   stress and represents the square root of the **elastic strain energy
%   density** per unit volume, restricted to the positive energy domain.
%
%--------------------------------------------------------------------------
% INPUTS:
%   STRAIN   : Strain tensor(s) in Voigt notation, possibly stacked for
%              several Gauss points (size = nstrain × nGAUSS, vectorized as
%              [ε₁; ε₂; ...]).
%
%   STRESSeff: Effective stress tensor(s) in Voigt notation, same layout as
%              STRAIN (size = nstrain × nGAUSS, vectorized).
%
%   nstrain  : Number of strain/stress components per Gauss point:
%                 • 3 → plane stress / 1D bar
%                 • 4 → plane strain (Voigt order: [xx, yy, xy, zz])
%                 • 6 → full 3D case
%
%--------------------------------------------------------------------------
% OUTPUTS:
%   r        : (vector, size nGAUSS×1)
%              Strain-like internal variable driving the damage law at each
%              Gauss point.
%
%--------------------------------------------------------------------------
% ALGORITHM:
%   1. Loop over strain/stress components.
%   2. Compute the scalar product ε : σ_eff = Σ_i ε_i σ_eff,i.
%   3. Take the square root to obtain the energy-equivalent measure r.
%
%--------------------------------------------------------------------------
% REMARKS:
%   • The formulation assumes that both STRAIN and STRESSeff are arranged
%     in contiguous blocks corresponding to each Gauss point.
%   • This definition preserves symmetry between strain and stress and is
%     suitable for scalar isotropic damage models (see Simo & Ju, 1987).
%   • The function is used in `DamageSmallStrainLargeRotations` to determine
%     whether each Gauss point remains elastic or undergoes damage evolution.
%
%--------------------------------------------------------------------------
% REFERENCES:
%   Simo, J.C. & Ju, J.W., "Strain- and stress-based continuum damage models:
%   I. Formulation and II. Computational aspects", Int. J. Solids Struct., 1987.
%
%--------------------------------------------------------------------------
% HISTORY:
%   - 21-Oct-2025, J.A. Hernández Ortega
%     First implementation for energy-equivalent r variable.
%--------------------------------------------------------------------------
%==========================================================================



nGAUSS= size(STRAIN,1)/nstrain ; 
r = zeros(nGAUSS,1) ; 

for istrain = 1:nstrain 
    INDloc = istrain:nstrain:size(STRAIN,1);
    r = r + STRAIN(INDloc).*STRESSeff(INDloc) ; 
end

r = sqrt(r) ; 