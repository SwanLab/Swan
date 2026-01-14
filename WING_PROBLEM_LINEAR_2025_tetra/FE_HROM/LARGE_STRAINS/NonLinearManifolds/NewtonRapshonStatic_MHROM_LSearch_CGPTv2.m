function [VAR,CONVERGED,DATA] = NewtonRapshonStatic_MHROM_LSearch_CGPTv1(DATA,OPERFE,VAR,MATPRO,DOFl)
% NewtonRapshonStatic_MHROM_LSearch is a variation of NewtonRapshonStaticLarge_manifoldMDw
% It features line search
% JAHO, 29-Oct-2025, Wednesday, UPC,Terrassa

% NewtonRapshonStaticLarge_manifoldMDw is a variation of
% NewtonRapshonStaticLarge_manifoldMD, described below
% The goal of the variation is to account for integration weights which
% varies with latent variables
% JAHO, 22-Sept-2025, HGs Pedralbes, Barcelona
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/11_MAW_ECM_plast.mlx
%--------------------------------------------------------------------------
% NewtonRapshonStaticLarge_manifoldMD
%--------------------------------------------------------------------------
% Solve static equilibrium on a nonlinear reduced manifold u = τ(q) using a
% Newton–Raphson scheme posed in the manifold’s tangent space. This MD
% variant supports **multi-latent** reduced coordinates (vector q).
%
% SYNTAX
%   [VAR, CONVERGED, DATA] = NewtonRapshonStaticLarge_manifoldMD( ...
%       DATA, OPERFE, VAR, MATPRO, DOFl)
%
% DESCRIPTION
%   At each Newton iteration:
%     1) Map reduced coords to full space with τ(q): u = [τ(qL); qR].
%     2) Assemble full residual/internal forces from u.
%     3) Project (pull back) residual and external load to reduced space via τ′(q):
%            r_ROM = τ′(q)^T r(u),   f_ROM = τ′(q)^T f_ext(u).
%     4) Assemble full tangent K(u) and project:
%            K_ROM = τ′(q)^T K(u) τ′(q).
%     5) (Optional) Add geometric term from τ″(q)·r(u) if τ″ available.
%     6) Solve K_ROM(DOFl,DOFl) ΔqL = −r_ROM(DOFl) and update qL.
%     7) Reconstruct fields, update energies/stresses, and test convergence.
%
% INPUTS
%   DATA   : Struct with solver controls and offline τ-evaluator.
%            .NEWTON_RAPHSON.NMAXiter   - Max Newton iterations
%            .ISDYNAMIC (0/1)           - If 1, dynamic residual/Jacobian (not adapted here)
%            .SMALL_STRAIN_KINEMATICS   - 0: large strain, 1: small strain
%            .DATA_evaluateTAU_and_DER  - Precompiled τ/τ′/τ″ evaluator descriptor
%            .ListFieldInternalVariables- Names of internal variables to carry over
%            .INTERNAL_FORCES_USING_precomputed_Kstiff - Flag for internal force path
%            .kiter                     - (set internally) current Newton iter
%            (… plus standard energy/tolerance fields used by CheckConvergenceLSTR)
%
%   OPERFE : FE operators/partitions and optional extras.
%            .DOFr, .DOFm, A, G, uBAR   - Partitioning/maps for constrained DOFs
%            .HYDRO (optional)          - Hydrostatic coupling operators
%            .wSTs_cluster (optional)   - Clustered ECM weight schedule with:
%                 .IndexDOFl_q          - Indices of qL within VAR.DISP
%                 .q                     - Grid of latent abscissae
%                 .Values                - Weight vectors per abscissa
%
%   VAR    : State container (updated in place).
%            .DISP                      - Reduced coords (qL on DOFl; qR on DOFr)
%            .FEXT_extended             - External forces in extended space
%            .(internal variables)      - As listed in DATA.ListFieldInternalVariables
%            (Fields updated inside: FINT, RESID, stresses, energies, etc.)
%
%   MATPRO : Constitutive parameters passed to residual/tangent assembly.
%
%   DOFl   : Indices of free reduced DOFs (subset of VAR.DISP corresponding to qL).
%
% OUTPUTS
%   VAR       : Updated state upon exit. On convergence includes Cauchy/Von Mises
%               stresses, energies, and last residuals/forces in reduced space.
%   CONVERGED : 1 if convergence criteria satisfied; 0 otherwise.
%   DATA      : Updated with iteration snapshots/history (e.g., via UpdateIterativeSnapshot).
%
% KEY FEATURES (vs. FST variant)
%   • Multi-latent qL support with consistent τ/τ′/τ″ projections.
%   • Efficient, precompiled τ-evaluator to avoid anonymous-function overhead.
%   • Optional geometric stiffness contribution via τ″(q)·r(u) (if τ″ provided).
%   • ECM rule (weights) **auto-switching** near manifold abscissae when Newton stalls.
%
% CONVERGENCE & CHECKS
%   • Convergence evaluated on reduced DOFs using CheckConvergenceLSTR
%     with (FINT, FEXT, RESID) projected/pulled back to reduced space.
%   • If negative Jacobians (e.g., detF ≤ 0) are detected and the code path
%     disallows proceeding, the routine exits with CONVERGED = 0.
%   • Dynamic terms (mass/damping) are stubbed: DATA.ISDYNAMIC==1 triggers
%     an error (“not adapted yet”) to avoid silent misuse.
%
% GEOMETRIC TERM (optional)
%   If τ″(q) is available, the projected tangent for DOFl is augmented by
%       K_geo = Σ_i ( τ″_i(q) * r_i(u) ),
%   where i runs over generalized components in the extended displacement
%   space. This captures curvature of the manifold in the Newton update.
%
% BOUNDARY CONDITIONS
%   • Unconstrained case: direct solve on K(DOFl,DOFl).
%   • With constraints: use condensed system via A, and then reconstruct
%     constrained DOFs by u_r = G u_m + ū.
%
% ECM WEIGHT SCHEDULING (when OPERFE.wSTs_cluster is provided)
%   On reaching the max Newton iterations, the routine can switch to the
%   nearest set of clustered ECM weights based on the current qL location,
%   reset the Newton counter, and retry (up to a small max number of retries).
%   This mitigates localization mismatches across manifold patches.
%
% PERFORMANCE TIPS
%   • Ensure τ/τ′/τ″ evaluators are vectorized and return sparse matrices
%     when appropriate (τ′ often is tall-skinny but structured).
%   • Use iterative solvers (PCG) only with good preconditioners; otherwise
%     direct solves are usually more robust in ROMs of modest size.
%   • If eigenvalue checks are enabled for diagnostics, guard them with flags
%     to avoid overhead in production runs.
%
% ASSUMPTIONS / LIMITATIONS
%   • Static formulation (dynamic contributions intentionally disabled).
%   • Large-strain support relies on ResidualFromDisplacementsVAR and
%     GetStiffnessMatrix providing consistent PK1/Cauchy and geometric terms.
%   • τ″(q) optional; setting it empty skips geometric correction.
%   • Hydrostatic forces path requires OPERFE.HYDRO to be fully defined.
%
% POSSIBLE FAILURE MODES
%   • Singular/ill-conditioned K(DOFl,DOFl): manifold region with weak tangent
%     or physical buckling; consider line search, damping, or ECM re-weighting.
%   • Persistent residual stagnation: try enabling τ″(q), refine τ offline, or
%     broaden ECM weights (cluster switch already attempted here).
%
% SEE ALSO
%   ResidualFromDisplacementsVAR, GetStiffnessMatrix, CheckConvergenceLSTR,
%   VonMisesCauchyStresses, UpdateVelocityAcceleration, ...
%   KineticAndStrainEnergyGet_mnfold
%
% REFERENCES (internal numbering)
%   MLEARNstruct_1.pdf, §9.2:
%     (176) τ(q) map; (182) residual projection; (184) tangent projection;
%     (185–193) geometric/dynamic extensions.
%
% HISTORY
%   • 02-Jul-2025  Original single-latent version (Barcelona).
%   • 11-Aug-2025  FST variant with fast τ evaluator (Cartagena).
%   • 13-Aug-2025  MD version (multi-latent q) and ECM auto-switch.
%   • 29-Aug-2025  Comment overhaul (ChatGPT).
%
% AUTHOR
%   Joaquín A. Hernández Ortega (UPC/CIMNE)
%--------------------------------------------------------------------------


if nargin == 0
    load('tmp.mat')
    %     CHECK_convergence =0;
    %     USE_GEOMETRIC_term_K = 0 ;
    %     USE_WEIGHT_TERM_k = 0 ;
    %USE_WEIGHT_TERM_k = 0 ;
end


% 
% kiter = 1;
% CONVERGED = 0 ;
% normd = 1e20 ;
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

% NR_MHROM_LineSearch_Reg
% -------------------------------------------------------------------------
% Robust Newton–Raphson with backtracking line search and progressive,
% eigen-aware regularization for a 2x2 latent system (MHROM-style).
%
% INPUTS
%   OPERFE, VAR, VARint_n, DATA, MATPRO, DOFl : your usual solver structs.
%   fun : struct with required function handles:
%       fun.GetResidual(OPERFE,VAR,VARint_n,DATA,MATPRO,DOFl,kiter)
%           -> [FintL,Rl,VAR,OPERFE,q_LAT,tauNON_q,tauNONder_q,tauNONder2_q, ...
%               celastST,FgradST,detFgrad,Kred_w_LL,DATA,VAR_FEXT_DOFl,tauNONallDER_q]
%       fun.GetJacobian(OPERFE,DATA,VAR,FgradST,celastST,tauNONallDER_q,tauNONder2_q,DOFl,Kred_w_LL)
%           -> Kll
%       UpdateVariablesNR(DATA,VAR,FgradST,detFgrad,OPERFE)
%           -> VAR
%       fun.CheckConvergence(FintL,VAR_FEXT_DOFl,Rl,kiter,normd,DATA,q_LAT)
%           -> [~, ~, CONVERGED]
%
% OUTPUTS
%   VAR  : state after solve (accepted last state).
%   info : struct with diagnostics:
%       .nIter, .finalResidual, .acceptedFullSteps, .nBacktracksTotal,
%       .muHistory, .lamMinHistory, .resHistory
%   DATA : returned in case you want to keep updated REG/LINESEARCH fields.
%
% Author: (you)
% -------------------------------------------------------------------------

% ---------- Defaults (safe) ----------
if ~isfield(DATA,'NEWTON_RAPHSON'), DATA.NEWTON_RAPHSON = struct; end
if ~isfield(DATA.NEWTON_RAPHSON,'NMAXiter'), DATA.NEWTON_RAPHSON.NMAXiter = 100; end
if ~isfield(DATA.NEWTON_RAPHSON,'TOL_RES'),  DATA.NEWTON_RAPHSON.TOL_RES  = 1e-6; end
if ~isfield(DATA.NEWTON_RAPHSON,'TOL_DU'),   DATA.NEWTON_RAPHSON.TOL_DU   = 1e-10; end
if ~isfield(DATA.NEWTON_RAPHSON,'PRINT'),    DATA.NEWTON_RAPHSON.PRINT    = true; end

if ~isfield(DATA,'LINESEARCH'), DATA.LINESEARCH = struct; end
if ~isfield(DATA.LINESEARCH,'c1'),             DATA.LINESEARCH.c1  = 1e-2; end   % Armijo
if ~isfield(DATA.LINESEARCH,'shrink'),         DATA.LINESEARCH.shrink = 0.5; end
if ~isfield(DATA.LINESEARCH,'max_backtracks'), DATA.LINESEARCH.max_backtracks = 10; end
if ~isfield(DATA.LINESEARCH,'min_step'),       DATA.LINESEARCH.min_step = 1e-8; end

if ~isfield(DATA,'REG'), DATA.REG = struct; end
if ~isfield(DATA.REG,'mu'),        DATA.REG.mu        = 0.0;  end % live LM parameter
if ~isfield(DATA.REG,'mu_max'),    DATA.REG.mu_max    = 1e4;  end
if ~isfield(DATA.REG,'pd_floor'),  DATA.REG.pd_floor  = 1e-8; end % target min eigen
if ~isfield(DATA.REG,'ramp_up'),   DATA.REG.ramp_up   = 2.0;  end % multiplicative grow
if ~isfield(DATA.REG,'ramp_dn'),   DATA.REG.ramp_dn   = 0.5;  end % multiplicative decay
if ~isfield(DATA.REG,'kappa'),     DATA.REG.kappa     = 1.0;  end % additive lift factor

% ---------- Diagnostics ----------
resHistory = zeros(DATA.NEWTON_RAPHSON.NMAXiter,1);
muHistory  = zeros(DATA.NEWTON_RAPHSON.NMAXiter,1);
lamHistory = zeros(DATA.NEWTON_RAPHSON.NMAXiter,1);

% ---------- Bookkeeping ----------
kiter = 1;
resf_old = inf;
prev_backtracks = 0;
last_step_accepted_full = false;

best_state.Rnorm = inf;

% ---------- Newton loop ----------
while kiter <= DATA.NEWTON_RAPHSON.NMAXiter
    
    % Residual etc. at CURRENT accepted state
    [FintL,Rl,VAR,OPERFE,q_LAT,~,tauNONder_q,tauNONder2_q,celastST,FgradST, ...
        detFgrad,Kred_w_LL,DATA,VAR_FEXT_DOFl,tauNONallDER_q] = ...
        GetResidual_MHROM(OPERFE,VAR,VARint_n,DATA,MATPRO,DOFl,kiter);
    
    if isempty(celastST) && isfield(DATA,'INTERNAL_FORCES_USING_precomputed_Kstiff') ...
            && DATA.INTERNAL_FORCES_USING_precomputed_Kstiff == 0
        if DATA.NEWTON_RAPHSON.PRINT
            disp('Not converged: Negative Jacobians in elements. Consider PROCEED_WITH_NEGATIVE_JACOBIANS=1.')
        end
        break
    end
    
    % Convergence check at current state
    normd_dummy = 1e20; % we haven’t computed step yet
    [~,~,CONVERGED] = CheckConvergenceLSTRq(FintL,VAR_FEXT_DOFl,Rl,kiter,normd_dummy,DATA,q_LAT);
    resf = norm(Rl);
    resHistory(kiter) = resf;
    
    if CONVERGED || (resf <= DATA.NEWTON_RAPHSON.TOL_RES)
        VAR = UpdateVariablesNR(DATA,VAR,FgradST,detFgrad,OPERFE);
        break
    end
    
    % Tangent
    Kll = GetJacobianMatrix_MHROM(OPERFE, DATA, VAR, FgradST, celastST, tauNONallDER_q, tauNONder2_q, DOFl, Kred_w_LL);
    
    % Eigen diagnostics (2x2, but robust for general small systems)
    Ksym = 0.5*(Kll + Kll.');
    lam = eig(Ksym);
    lam_min = min(lam);
    lamHistory(kiter) = lam_min;
    
    % Progressive regularization:
    % Additive lift if near singular; multiplicative ramp based on health/LS behavior
    healthy_curvature = lam_min > 10*DATA.REG.pd_floor;
    near_singular     = lam_min <= DATA.REG.pd_floor;
    
    mu_add = 0;
    if near_singular
        mu_add = DATA.REG.kappa * max(0, DATA.REG.pd_floor - lam_min);
    end
    
    mu_target = DATA.REG.mu; % start from current
    
    if near_singular
        mu_target = max(DATA.REG.mu + mu_add, DATA.REG.mu * DATA.REG.ramp_up);
    elseif (~last_step_accepted_full) || (prev_backtracks > 0)
        % line search struggled last iteration -> softly ramp up
        mu_target = max(DATA.REG.mu, DATA.REG.mu * DATA.REG.ramp_up);
    else
        % healthy curvature & easy full step last iter -> decay
        mu_target = DATA.REG.mu * DATA.REG.ramp_dn;
    end
    
    mu_target = min(mu_target, DATA.REG.mu_max);
    if mu_target < 1e-14, mu_target = 0; end
    DATA.REG.mu = mu_target;
    muHistory(kiter) = DATA.REG.mu;
    
    % Build regularized tangent
    Kreg = Kll + DATA.REG.mu * speye(size(Kll));
    
    % Search direction
    delta_dir = - Kreg \ Rl;
    normd = norm(delta_dir);
    
    % Backtracking line search (Armijo on ||R||)
    alpha = 1.0;
    R0_norm = resf;
    accepted = false;
    n_back = 0;
    
    % reset best-so-far for this iteration
    best_try = struct('Rnorm', inf);
    
    for ib = 1:DATA.LINESEARCH.max_backtracks
        VAR_try = VAR;
        VAR_try.DISP(DOFl) = VAR.DISP(DOFl) + alpha*delta_dir;
        
        % Enforce constraints if present
        if isfield(OPERFE,'DOFm') && ~isempty(OPERFE.DOFm)
            VAR_try.DISP(OPERFE.DOFr) = OPERFE.G * VAR_try.DISP(OPERFE.DOFm) + OPERFE.uBAR;
        end
        
        % Evaluate residual at trial state
        [~, Rl_try, VAR_try_loc, OPERFE_try, ~, ~, ~, ~, celastST_try, FgradST_try, ...
            detFgrad_try, ~, DATA_try, VAR_FEXT_DOFl_try, tauNONallDER_q_try] = ...
            GetResidual_MHROM(OPERFE, VAR_try, VARint_n, DATA, MATPRO, DOFl, kiter);
        
        % Reject if invalid trial (e.g., negative Jacobians and not allowed)
        if isempty(celastST_try) && isfield(DATA,'INTERNAL_FORCES_USING_precomputed_Kstiff') ...
                && DATA.INTERNAL_FORCES_USING_precomputed_Kstiff == 0
            alpha = alpha * DATA.LINESEARCH.shrink;
            n_back = n_back + 1;
            if alpha < DATA.LINESEARCH.min_step, break; end
            continue
        end
        
        R_try_norm = norm(Rl_try);
        
        % Track best-so-far
        if R_try_norm < best_try.Rnorm
            best_try = struct('Rnorm',R_try_norm,'alpha',alpha, ...
                'VAR',VAR_try_loc,'OPERFE',OPERFE_try,'FgradST',FgradST_try, ...
                'detFgrad',detFgrad_try,'celastST',celastST_try,'DATA',DATA_try, ...
                'VAR_FEXT_DOFl',VAR_FEXT_DOFl_try,'tauNONallDER_q',tauNONallDER_q_try);
        end
        
        % Armijo sufficient decrease: ||R(alpha)|| <= (1 - c1*alpha) * ||R(0)||
        if R_try_norm <= (1 - DATA.LINESEARCH.c1*alpha) * R0_norm
            % Accept
            VAR            = best_try.VAR;   % or VAR_try_loc (best_try == current typically)
            OPERFE         = best_try.OPERFE;
            FgradST        = best_try.FgradST;
            detFgrad       = best_try.detFgrad;
            celastST       = best_try.celastST;
            DATA           = best_try.DATA;
            VAR_FEXT_DOFl  = best_try.VAR_FEXT_DOFl;
            tauNONallDER_q = best_try.tauNONallDER_q;
            
            resf = best_try.Rnorm;
            accepted = true;
            last_step_accepted_full = (alpha >= 0.99);
            break
        else
            alpha = alpha * DATA.LINESEARCH.shrink;
            n_back = n_back + 1;
            if alpha < DATA.LINESEARCH.min_step
                break
            end
        end
    end
    
    if ~accepted
        % Fall back to best improvement if any
        if best_try.Rnorm < R0_norm
            VAR            = best_try.VAR;
            OPERFE         = best_try.OPERFE;
            FgradST        = best_try.FgradST;
            detFgrad       = best_try.detFgrad;
            celastST       = best_try.celastST;
            DATA           = best_try.DATA;
            VAR_FEXT_DOFl  = best_try.VAR_FEXT_DOFl;
            tauNONallDER_q = best_try.tauNONallDER_q;
            resf           = best_try.Rnorm;
        end
        last_step_accepted_full = false;
    end
    
    prev_backtracks = n_back;
    
    % Stagnation / safety exits
    if resf <= DATA.NEWTON_RAPHSON.TOL_RES
        VAR = UpdateVariablesNR(DATA,VAR,FgradST,detFgrad,OPERFE);
        break
    end
    if normd <= DATA.NEWTON_RAPHSON.TOL_DU
        if DATA.NEWTON_RAPHSON.PRINT
            disp('NR: Displacement update below TOL_DU. Stopping.');
        end
        VAR = UpdateVariablesNR(DATA,VAR,FgradST,detFgrad,OPERFE);
        break
    end
    if kiter > 5 && isfinite(resf_old) && abs(resf - resf_old)/max(1,resf_old) < 1e-4
        if DATA.NEWTON_RAPHSON.PRINT
            disp('NR: Stagnation detected. Stopping with best-so-far.');
        end
        VAR = UpdateVariablesNR(DATA,VAR,FgradST,detFgrad,OPERFE);
        break
    end
    
    % Print trace
    if DATA.NEWTON_RAPHSON.PRINT
        fprintf('iter=%3d  ||R||=%.3e  mu=%.3e  lam_min=%.3e  back=%d  full=%d\n', ...
            kiter, resf, DATA.REG.mu, lam_min, n_back, last_step_accepted_full);
    end
    
    resf_old = resf;
    kiter = kiter + 1;
end

% Finalize
info.nIter             = kiter;
info.finalResidual     = resfHistorySafe(resHistory, kiter);
info.acceptedFullSteps = last_step_accepted_full;
info.nBacktracksTotal  = sumBacktracks(prev_backtracks); % trivial here, keep for API
info.muHistory         = muHistory(1:kiter);
info.lamMinHistory     = lamHistory(1:kiter);
info.resHistory        = resHistory(1:kiter);

end % function

% ----------------- helpers -----------------
function r = resfHistorySafe(H, k)
if k<=0 || k>numel(H) || ~isfinite(H(k)), r = NaN; else, r = H(k); end
end

function s = sumBacktracks(b)
% Placeholder for future aggregation across iterations
s = b;
end
