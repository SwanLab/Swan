function  [Fint,fINTredNON_mst_DOFl,wMSTeq,ResREDeqMST,der_fINTredNON_slv_DOFl] = InternalForcesECMnon2md(OPERFE,PoneST,PK2STRESS,DATA,VAR,tauNONallDER_q,...
    PoneST_incre,PK2STRESS_incre) 
 % InternalForcesECMnon2md
%--------------------------------------------------------------------------
% Adaptation of InternalForc2sECMnon2 to cope with:
%   (i) two latent variables (typical of elastoplasticity: e.g., plastic strain
%       amplitude and hardening-like scalar), and
%   (ii) an additive split of stresses into linear (elastic predictor) and
%        nonlinear (plastic correction) parts.
%
% JAHO, 29 Aug 2025, Balmes 185 (Barcelona)
%
% PURPOSE
% -------
% Compute the reduced internal force vector Fint for a nonlinear HROM on a
% τ(q)-manifold where the internal force density is assembled from:
%   • MASTER (ECM) contribution: evaluated at selected “master” quadrature
%     points using first Piola–Kirchhoff stresses (P1 = PoneST).
%   • SLAVE contribution: approximated/extrapolated from the master part via
%     a nonlinear surrogate map η_NON that depends on the manifold tangent
%     τ′(q) and *two* latent variables (e.g., s₁, s₂) tracking the EP state.
%
% The routine also returns equivalent master weights and residuals used in
% tangent linearization and Newton steps on the manifold.
%
% CONTEXT / MODELING NOTES
% -----------------------
% • Manifold mechanics: State variables live on a reduced manifold τ(q). The
%   mapping τ and its Jacobian τ′(q) are precompiled; τ′(q) is provided here
%   as tauNONallDER_q to project internal work densities onto the tangent space.
%
% • Two latent variables (elastoplasticity):
%   The surrogate η_NON(·) internally carries two latent coordinates that
%   distinguish elastic/plastic branches and loading/unloading sub-branches.
%   This enables *piecewise-injective* behavior with smooth transitions across
%   branches, avoiding non-injectivity of master–slave mappings in localized EP.
%
% • Linear/nonlinear stress split:
%     P = P^lin + P^non,   with   P^lin = C₀ : E  (elastic predictor),
%   while P^non captures plastic/algorithmic corrections (return-mapping).
%   The master evaluation uses the supplied PoneST (P), and the slave part is
%   reconstructed by η_NON from the master internal force densities projected
%   with τ′(q). This keeps linear terms exact where desired and confines
%   approximation to the genuinely nonlinear share.
%
% • Fallback (PK2 path):
%   If PoneST is empty, the function offers a special path for PK2 stresses
%   (mainly for EIFEM/testing). The general non-CECM path is *not yet*
%   implemented for this md-variant (two-latent + split).
%
% INPUTS
% ------
% OPERFE : struct with offline/assembly operators for ECM/HROM:
%   .Bst               Stress projection operator at ECM points (master DOFs)
%   .wSTs              Weights at master ECM points (vector matching PoneST blocks)
%   .DOFl              Free-DOF indices in reduced/global layout
%   .DATA_regress_eta_der  (if present) precompiled evaluator for η_NON and its derivatives
%   .wRED_slv          Weights/coefficients used by the slave extrapolation
%   .KstiffLINEAR      (optional) linear stiffness (used in special PK2 path)
%
% PoneST : [ngaus_master*ndim^2 × 1] or [] 
%   Stacked 1st Piola–Kirchhoff stresses at master ECM points. When provided,
%   the “LT” path is used: InternalForcesECMnon2LT(...).
%
% PK2STRESS : [ngaus_master*nstrain × 1] or [] 
%   Stacked 2nd Piola–Kirchhoff stresses for the fallback “SP” path used when
%   PoneST is empty and DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES==1.
%
% DATA : struct
%   .MESH.ndim         Spatial dimension (2 or 3)
%   .MESH.nstrain      Number of strain components
%   .CECM_ONLY_FOR_NONLINEAR_STRESSES  Flag controlling the PK2 fallback mode
%
% VAR : struct
%   Current step state (e.g., external forces, optional incremental PK2 for
%   EIFEM special path, etc.). Exact fields depend on your codebase.
%
% tauNONallDER_q : [n_modes_all × n_modes_manifold]
%   Jacobian τ′(q) of the manifold map at the current (reduced) configuration.
%   Used to project master internal work densities onto the manifold tangent.
%
% OUTPUTS
% -------
% Fint : [n_modes_all × 1]
%   Reduced internal force vector combining master ECM evaluation and
%   slave-domain extrapolation through η_NON (two-latent, split-aware).
%
% fINTredNON_mst_DOFl : [n_modes_manifold × n_master_points]
%   Reduced internal *force density* at master points, projected to the
%   manifold tangent and restricted to free DOFs (DOFl). This is the key
%   quantity fed to η_NON to reconstruct the slave contribution.
%
% wMSTeq : [n_master_points × 1]
%   Equivalent (nonlinear) master weights incorporating the extrapolated
%   slave share; they enable expressing the full internal force as a
%   master-only quadrature (cf. theory Eq. (9.19)/(9.29)-style expressions).
%
% ResREDeqMST : [n_modes_all × 1]
%   Residual based on the equivalent master integration. Useful to assemble
%   tangents/Jacobians consistently with the η_NON-enhanced quadrature
%   (cf. theory Eq. (9.33)-type residual).
%
% der_fINTredNON_slv_DOFl : struct / array
%   Derivatives needed for tangent linearization of the *slave* contribution:
%     • w.r.t. manifold DOFs via τ′(q) (and, internally, latent variables),
%     • optionally w.r.t. η_NON parameters if using adaptive/learned maps.
%   Exact shape/fields follow the implementation in InternalForcesECMnon2LT/SP.
%
% CONTROL FLOW
% ------------
% if ~isempty(PoneST)
%   → Call InternalForcesECMnon2LT(...)  % “LT” = linear/nonlinear split aware,
%                                        % two-latent η_NON path using P1.
% else
%   if DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
%       error('Option not implemented yet, 29-August-2025')
%   else
%       → Call InternalForcesECMnon2SP(...)  % Special PK2/EIFEM-oriented path
%                                            % retaining the two-latent handling.
%   end
% end
%
% ASSUMPTIONS & REQUIREMENTS
% -------------------------
% • PoneST must be block-ordered consistently with OPERFE.wSTs and Bst.
% • tauNONallDER_q corresponds to the *current* q and matches the reduced
%   bases used to build OPERFE operators (dimensionally compatible).
% • η_NON and its derivatives (inside the LT/SP workers) are trained/fitted
%   with two-latent inputs and support the linear/nonlinear stress split.
%
% NUMERICAL REMARKS
% -----------------
% • The two-latent, branch-aware η_NON improves robustness in localized EP
%   (plastic bands) where injectivity of a single master–slave map fails.
% • The split avoids polluting elastic predictor contributions with
%   approximation error; only plastic (nonlinear) corrections are extrapolated.
% • For tangent consistency, always use the returned der_fINTredNON_slv_DOFl
%   together with wMSTeq/ResREDeqMST in your Newton assembly.
%
% FAILURE MODES / DIAGNOSTICS
% ---------------------------
% • error('Option not implemented...'): raised when attempting the non-CECM
%   PK2 path with DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES==0 in the md variant.
% • Dimensional mismatches typically originate from inconsistent stacking of
%   stress vectors (P1 vs. PK2) or τ′(q) dimension errors—check sizes early.
%
% VERSION HISTORY
% ---------------
% • 2025-08-29 (JAHO): first md-variant with two latent vars + split handling.
% • 2025-07-22 (JAHO): base InternalForcesECMnon2 documented (single latent).
%
% GLOSSARY (EN → ES)
% ------------------
% • latent variables → variables latentes
% • master/slave domain → dominio maestro/esclavo
% • split (linear/nonlinear) → descomposición (lineal/no lineal)
% • equivalent weights → pesos equivalentes
% • tangent (linearization) → tangente (linealización)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

%
if nargin == 0
    load('tmp2.mat')
end


if ~isempty(PoneST)
    [Fint,fINTredNON_mst_DOFl,wMSTeq,ResREDeqMST,der_fINTredNON_slv_DOFl] = InternalForcesECMnon2LT(OPERFE,PoneST,DATA,VAR,tauNONallDER_q)  ;
else
    
    if  DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
        error('Option not implemented yet, 29-August-2025')
        
        nF = DATA.MESH.nstrain ;
        for icomp = 1:nF
            icol = icomp:nF:length(PK2STRESS) ;
            PK2STRESS(icol,:) = PK2STRESS(icol,:).*OPERFE.wSTs;
        end
        Fint = OPERFE.Bst'*PK2STRESS ;
    else
        %         % Special implementation for EIFEM
        %         % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
        %         nF = DATA.MESH.nstrain ;
        %         for icomp = 1:nF
        %             icol = icomp:nF:length(PK2STRESS) ;
        %             PK2STRESS(icol,:) = VAR.PK2STRESS_incre(icol,:).*OPERFE.wSTs;
        %         end
        %         Fint = OPERFE.KstiffLINEAR*VAR.DISP + OPERFE.Bst'*PK2STRESS ;
        
        [Fint,fINTredNON_mst_DOFl,wMSTeq,ResREDeqMST,der_fINTredNON_slv_DOFl] = InternalForcesECMnon2SP(OPERFE,PK2STRESS_incre,DATA,VAR,tauNONallDER_q)  ;
        
        
    end
end

