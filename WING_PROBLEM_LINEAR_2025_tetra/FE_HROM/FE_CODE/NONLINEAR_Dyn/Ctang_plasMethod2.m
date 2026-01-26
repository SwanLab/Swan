function CtangPlas  =  Ctang_plasMethod2(nstrain,normSTRgA,pmultA,normalNP,muG,kappaG,HmodulusG,identVOL,identV) 
% =========================================================================
% Ctang_plasMethod2 — Consistent tangent for small-strain J2 plasticity (vectorized)
% =========================================================================
% PURPOSE
%   Build the *algorithmic (consistent) tangent* C_t for a small-strain,
%   pressure-dependent split (volumetric + deviatoric) formulation of
%   associative J2 plasticity with **isotropic hardening**. The routine is
%   fully vectorized over Gauss points and returns the block-stacked
%   Voigt form used downstream to assemble the global tangent.
%
% FORM OF THE TANGENT (per Gauss point)
%     C_t = λ_pl I⊗I + 2 μ_pl I_dev + ρ_N (n ⊗ n)
%
%   where
%     - λ_pl, μ_pl : plastic-corrected Lamé moduli
%     - I⊗I       : volumetric projector (volumetric part)
%     - I_dev     : deviatoric projector in Voigt (with 0.5 on γ_xy row/col)
%     - n         : outward normal to the yield surface in Voigt space
%     - ρ_N       : scalar weight multiplying the dyadic n⊗n (plastic part)
%
% INPUTS
%   nstrain     : number of Voigt components (e.g., 4 for plane strain: [xx yy xy zz])
%   normSTRgA   : ||s_tr||, norm of the deviatoric *trial* stress (vector, per GP)
%   pmultA      : plastic multiplier increment Δγ (vector, per GP)
%   normalNP    : yield normal n in Voigt form, stacked for all GPs
%                 (length = nstrain * ngaus). Layout: [n_1; n_2; ...; n_ngaus]
%                 with each n_i of length nstrain.
%   muG         : shear modulus μ at each GP (vector, length = ngaus)
%   kappaG      : bulk modulus κ at each GP (vector, length = ngaus)
%   HmodulusG   : isotropic hardening modulus H at each GP (vector, length = ngaus)
%   identVOL    : block-stacked volumetric projector ( (ngaus*nstrain) × nstrain )
%   identV      : block-stacked deviatoric projector ( (ngaus*nstrain) × nstrain )
%
% OUTPUTS
%   CtangPlas   : block-stacked consistent tangent in Voigt form
%                 size = (ngaus*nstrain) × nstrain.
%
% INTERNAL DEFINITIONS (per GP)
%   betaPLAS   = 2 Δγ ( μ / ||s_tr|| )
%   μ_pl       = μ (1 − betaPLAS)
%   λ_pl       = κ − (2/3) μ_pl
%   coeffH     = 1 + H / (3 μ)
%   ρ_N        = 2 μ ( betaPLAS − 1/coeffH )
%
%   → Decomposition:
%       C_vol    = λ_pl I⊗I
%       C_dev    = 2 μ_pl I_dev
%       C_normal = ρ_N (n ⊗ n)
%
% VECTORIZATION / SHAPES
%   - identVOL and identV repeat the 4×4 (or nstrain×nstrain) projectors for
%     each Gauss point; bsxfun/@times scales rows by the per-GP scalars.
%   - C_normal is built column-by-column as ρ_N n n^T:
%       * ρ_N_N = ρ_N .* n  (Hadamard)
%       * For each Voigt component j, C_normal(:, j) = ρ_N_N .* n_j
%     where n is reshaped to (nstrain × ngaus), then re-stacked.
%
% ASSUMPTIONS / NOTES
%   - Small strain, associative J2 with linear isotropic hardening.
%   - normalNP must be the *unit* normal in the chosen Voigt metric
%     consistent with identV (i.e., engineering shear with the 0.5 factor).
%   - normSTRgA > 0 on plastic points; if very small, guard upstream to avoid
%     division by ~0 in betaPLAS.
%   - In elasticity (Δγ = 0) → betaPLAS = 0, ρ_N = −2μ(1 − 1/coeffH) ≤ 0;
%     however, this routine is intended for plastic GPs; pure elastic GPs
%     should use C = λ I⊗I + 2μ I_dev.
%
% UNITS
%   μ, κ, H in stress units (e.g., MPa). CtangPlas returns the same units.
%
% STABILITY TIPS
%   - If ||s_tr|| is tiny, cap betaPLAS using a small floor on normSTRgA.
%   - Ensure n is computed from the *trial* deviatoric stress direction and
%     normalized in the same inner product used by the return mapping.
%   - For robustness in strong softening (small H), consider line search or
%     consistent algorithmic updates of Δγ and n.
%
% REFERENCES
%   - Simo & Hughes, *Computational Inelasticity*, Springer.
%   - de Souza Neto et al., *Computational Methods for Plasticity*, Wiley.
% =========================================================================
% JAHO, commented on Oct 21st 2025, automatically by ChatGPT

if nargin == 0
    load('tmp3.mat')
end
% Coefficientes
betaPLAS =  2*pmultA.*(muG./normSTRgA);
muPLAS = muG.*(1-betaPLAS) ; 
lambdaPLAS = kappaG-2*muPLAS/3 ;
coeffH = 1+ HmodulusG./(3*muG) ; 
rhoN = 2*muG.*(betaPLAS-1./coeffH) ; 


%%%%
% Volumetric term 
C_vol =  bsxfun(@times,identVOL,lambdaPLAS) ; 
% Deviatoric term 
C_dev =  bsxfun(@times,identV,2*muPLAS) ; 
% Normal term 
C_normal= zeros(length(muG),nstrain); 
rhoN_N =  rhoN.*normalNP ;  
ngausG = length(normalNP)/nstrain ; 
normalNP_comp = reshape(normalNP,nstrain,ngausG) ; 
for istrain = 1:nstrain 
    nCOMPaux =  repmat(normalNP_comp(istrain,:),nstrain,1) ; 
    nCOMP =  reshape(nCOMPaux,nstrain*ngausG,1) ; 
    C_normal(:,istrain) = rhoN_N.*nCOMP ;
end

% Total 
CtangPlas = C_vol + C_dev + C_normal ;

 


% % Volumetric part
% KonePoneB = VolPartCtangBElas(BBfNW(indPLASm,:),nstrain,lambdaPLAS);
% %   part multiplying identity matrix 
% ngausT = length(indPLASm) ;
% muPLASa =reshape([muPLAS'; muPLAS'; muPLAS'; muPLAS'],ngausT,1) ; % For all gauss points and components
% dMuIdentB = CalcIdentB(BBfNW(indPLASm,:),nstrain,2*muPLASa)  ; 
% %   part multiplying nxn
% Npart = NPartCtang(BBfNW(indPLASm,:),nstrain,rhoN,normalNP);
% 
% 
% % Total 
% CtangB  = KonePoneB + dMuIdentB + Npart ;