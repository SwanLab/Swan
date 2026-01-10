function [VAR,CONVERGED,DATA] = NewtonRapshonStaticLarge_manifoldMDw(DATA,OPERFE,VAR,MATPRO,DOFl,USE_GEOMETRIC_term_K,USE_WEIGHT_TERM_k)
 %--------------------------------------------------------------------------
% NewtonRapshonStaticLarge_manifoldMDw
% -----------------------------------
% PURPOSE
%   Static Newton–Raphson solver on a nonlinear τ(q)-manifold (large-strain)
%   supporting **multi-latent** reduced coordinates with two optional
%   contributions:
%     • GEOMETRIC term in the reduced tangent (curvature of τ via τ″·r),
%     • WEIGHT term that accounts for **latent-dependent** ECM integration
%       weights (MAW-ECM setting).
%   This “MDw” variant is meant for cases where the hyperreduction rule
%   (sampling/weights) varies across the manifold abscissae, and you want to
%   explicitly include that dependence in the Jacobian assembly.
%
% SYNTAX
%   [VAR, CONVERGED, DATA] = NewtonRapshonStaticLarge_manifoldMDw( ...
%       DATA, OPERFE, VAR, MATPRO, DOFl, USE_GEOMETRIC_term_K, USE_WEIGHT_TERM_k)
%
% DESCRIPTION (per Newton iteration)
%   1) Residual assembly on the manifold:
%        [Fint_L, R_L, ..., celastST, FgradST, detFgrad, Kred_w_LL, ...]
%        = GetResidual_MHROM(...)
%      where Kred_w_LL is the **ECM-weighted** reduced tangent block built at
%      the current iterate (may already reflect latent-dependent weights).
%   2) Convergence check on DOFl:
%        CheckConvergenceLSTRq(Fint_L, Fext_L, R_L, kiter, normΔ, DATA, q_LAT)
%   3) Reduced Jacobian (with optional terms):
%        [Kll, Kll_geo] = GetJacobianMatrix_MHROM(..., DOFl, Kred_w_LL, ...
%                          USE_GEOMETRIC_term_K, USE_WEIGHT_TERM_k)
%      • If USE_GEOMETRIC_term_K = 1, include geometric (τ″·r) contribution.
%      • If USE_WEIGHT_TERM_k   = 1, include derivatives induced by changes
%        in ECM weights with respect to latent coordinates (MAW-ECM).
%   4) Newton update (no line search in this variant):
%        δqL = −Kll \ R_L
%        qL ← qL + δqL
%        If affine/periodic BCs: u_r ← G u_m + ū
%   5) On convergence: UpdateVariablesNR to finalize stresses/energies.
%
% INPUTS (summary)
%   DATA   : Solver controls (NEWTON_RAPHSON.*, tolerances/flags, negative-J handling, etc.).
%   OPERFE : Operators/partitions/maps (DOFr/DOFm, A, G, uBAR), ECM data.
%   VAR    : State container (DISP, FEXT_extended, internal vars listed in DATA.ListFieldInternalVariables).
%   MATPRO : Constitutive parameters.
%   DOFl   : Indices of free latent coordinates.
%   USE_GEOMETRIC_term_K : logical flag to include geometric tangent term.
%   USE_WEIGHT_TERM_k    : logical flag to include latent-dependent weight term.
%
% OUTPUTS
%   VAR       : Updated state (on success, consistent stresses/energies).
%   CONVERGED : 1 if convergence within NMAXiter; 0 otherwise.
%   DATA      : Iteration history/snapshots updated for diagnostics.
%
% KEY FEATURES / DIFFERENCES (vs. MD and GROK)
%   • MDw = MD + explicit switches for:
%       - Geometric stiffness contribution (τ″·r) in Kll.
%       - Latent-dependent ECM weight term in Kll (MAW-ECM sensitivity).
%   • No line search here (contrast with *GROK* variant which adds Armijo LS
%     plus spectral regularization).
%
% CONVERGENCE & SAFETY
%   • Uses CheckConvergenceLSTRq on DOFl.
%   • If negative Jacobians are detected and DATA disallows proceeding, exit
%     non-converged after snapshotting the iterate.
%   • Static solver; dynamic contributions are not considered here.
%
% PRACTICAL TIPS
%   • Keep GetResidual_MHROM and GetJacobianMatrix_MHROM consistent in
%     how ECM weights (and their latent sensitivity) are applied.
%   • If convergence is fragile, try enabling USE_GEOMETRIC_term_K first;
%     if residuals stall near patch boundaries, consider USE_WEIGHT_TERM_k
%     or switching to the GROK (line-search) variant.
%
% SEE ALSO
%   GetResidual_MHROM, GetJacobianMatrix_MHROM, CheckConvergenceLSTRq,
%   UpdateVariablesNR, NewtonRapshonStatic_MHROM_LSearch_GROK.
%
% REFERENCES (internal)
%   MLEARNstruct_1.pdf, §9.2: τ(q), projection; §9.3: ECM hyperreduction.
%   MAW-ECM notes: 112_NonLIN_ROM_RBF/11_MAW_ECM_plast.mlx
%
% HISTORY
%   • 02-Jul-2025  Original single-latent version (Barcelona).
%   • 11-Aug-2025  FST variant with fast τ evaluator (Cartagena).
%   • 13-Aug-2025  MD version (multi-latent q) and ECM auto-switch.
%   • 29-Aug-2025  Comment overhaul (ChatGPT).
%   • 22-Sept-2025 MDw variant for latent-dependent weights (HGs Pedralbes, Barcelona).
%
% AUTHOR
%   Joaquín A. Hernández Ortega (JAHO)
%
% COMMENTS UPDATED BY: ChatGPT-5
%   Date of update: 07-NOV-2025, Barcelona, Spain
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------


if nargin == 0
    load('tmp.mat')
    %     CHECK_convergence =0;
    %     USE_GEOMETRIC_term_K = 0 ;
    %     USE_WEIGHT_TERM_k = 0 ;
    %USE_WEIGHT_TERM_k = 0 ;
end



kiter = 1;
CONVERGED = 0 ;
normd = 1e20 ;
% Internal variables (value at previous time step)....Why do we need that?
% 
% ----------------------
VARint_n = [];
if ~isempty(DATA.ListFieldInternalVariables)
    for iINTVAR = 1:length(DATA.ListFieldInternalVariables)
        NameIntVar= DATA.ListFieldInternalVariables{iINTVAR} ;
        VARint_n.(NameIntVar) = VAR.(NameIntVar) ;
    end
end



while  kiter <=DATA.NEWTON_RAPHSON.NMAXiter
    
    [FintL,Rl,VAR,OPERFE,q_LAT,tauNON_q,tauNONder_q,tauNONder2_q,celastST,FgradST...
        ,detFgrad,Kred_w_LL,DATA,VAR_FEXT_DOFl,tauNONallDER_q]...
        = GetResidual_MHROM(OPERFE,VAR,VARint_n,DATA,MATPRO,DOFl,kiter) ;    
    
    if isempty(celastST) && DATA.INTERNAL_FORCES_USING_precomputed_Kstiff == 0
        disp('Not converged ..., because of Negative JAcobians (you may disable this by setting PROCEED_WITH_NEGATIVE_JACOBIANS=1)')
        DATA = ReshapeIterativeSnapshot(DATA,VAR,kiter) ;
        break
    else
        % 7. Check convergence        
        [~, ~,CONVERGED]=  CheckConvergenceLSTRq(FintL,VAR_FEXT_DOFl,Rl,kiter,normd,DATA,q_LAT) ;        
        
        if CONVERGED == 1
           VAR = UpdateVariablesNR(DATA,VAR,FgradST,detFgrad,OPERFE)  ; 
            break
        else            
            [Kll,Kll_geo] = GetJacobianMatrix_MHROM(OPERFE, DATA, VAR, FgradST, celastST,tauNONallDER_q,tauNONder2_q,DOFl,Kred_w_LL,...
                USE_GEOMETRIC_term_K,USE_WEIGHT_TERM_k)  ;
            
%             norm_Kll_mat = norm(Kll-Kll_geo-Kred_w_LL,'fro') ;
%             norm_Kll_w = norm(Kred_w_LL,'fro') ;
%              norm_Kll_geo = norm(Kll_geo,'fro') ;
%              norm_Kll = norm(Kll,'fro') ; 
                      
            % Search direction 
            delta_dL = -Kll\Rl;
            normd = norm(delta_dL) ;
            
            % Update displacement 
            VAR.DISP(DOFl) =  VAR.DISP(DOFl) + delta_dL ;
            
            if  ~isempty(OPERFE.DOFm)
                VAR.DISP(OPERFE.DOFr) =  OPERFE.G*VAR.DISP(OPERFE.DOFm) + OPERFE.uBAR ;
            end
            
            
        end
        kiter = kiter + 1;
        
        
        
        
    end
end