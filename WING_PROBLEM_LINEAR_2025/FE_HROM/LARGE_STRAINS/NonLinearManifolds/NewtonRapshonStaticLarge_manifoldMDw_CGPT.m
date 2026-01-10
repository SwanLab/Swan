function [VAR,CONVERGED,DATA] = NewtonRapshonStaticLarge_manifoldMDw_CGPT(DATA,OPERFE,VAR,MATPRO,DOFl,USE_GEOMETRIC_term_K,...
    USE_WEIGHT_TERM_k)
%NEWTONRAPSHONSTATICLARGE_MANIFOLDMDW  Robust NR solver on a nonlinear manifold u = τ(q).
% Rewritten with globalization and damping (LM + Armijo backtracking).
%
% SYNTAX
%   [VAR, CONVERGED, DATA] = NewtonRapshonStaticLarge_manifoldMDw( ...
%       DATA, OPERFE, VAR, MATPRO, DOFl, USE_GEOMETRIC_term_K, USE_WEIGHT_TERM_k)
%
% MAJOR CHANGES (vs. prior version)
%   • Adds LM regularization on reduced tangent K(DOFl,DOFl).
%   • Adds Armijo backtracking on merit φ = 1/2 ||r_DOFl||^2.
%   • Adds strong stopping tests (residual + step).
%   • Encapsulates trial evaluations to avoid code duplication.
%   • Keeps optional manifold-geometry and ECM-weight contributions.
%
% EXPECTED EXTERNALS (unchanged)
%   ResidualFromDisplacementsVARw, GetStiffnessMatrix,
%   VonMisesCauchyStresses, KineticAndStrainEnergyGet_mnfold,
%   UpdateIterativeSnapshot
%
% NOTES
%   • DATA.DATA_evaluateTAU_and_DER.nameFunctionEvaluate(qL, DESCR)
%       must return [tauNON_q, tauNONder_q, tauNONder2_q].
%   • OPERFE.DOFr/DOFm/A/G/uBAR are used when constraints are present.
%   • Hydro and dynamic paths are intentionally kept as in your base code.
%
% AUTHOR
%   Rewritten for robustness by ChatGPT for Joaquín A. Hernández Ortega (UPC/CIMNE)
%   28-Oct-2025, Barcelona

% -------------------------------------------------------------------------
% 0) Defaults & quick init
% -------------------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
    %     USE_GEOMETRIC_term_K = 0;
    %     USE_WEIGHT_TERM_k    = 0;
end

if ~isfield(DATA,'NEWTON_RAPHSON') || ~isfield(DATA.NEWTON_RAPHSON,'NMAXiter')
    DATA.NEWTON_RAPHSON.NMAXiter = 30;
end

% Globalization / damping knobs (safe defaults)
DATA = set_default(DATA,'LM_lambda',1e-6);        % initial LM damping
DATA = set_default(DATA,'LM_lambda_min',1e-14);
DATA = set_default(DATA,'LM_lambda_max',1e+8);
DATA = set_default(DATA,'armijo_c',1e-4);         % sufficient decrease
DATA = set_default(DATA,'armijo_rho',0.5);        % backtracking factor
DATA = set_default(DATA,'armijo_alpha_min',1e-6);

% Stopping tolerances (scaled)
DATA = set_default(DATA,'TOL_RES',1e-8);
DATA = set_default(DATA,'TOL_STEP',1e-10);

% Bookkeeping
kiter      = 1;
CONVERGED  = 0;
norm_step  = inf;

% Keep a copy of internal variables of the previous time step
VARint_n = pull_internal_variables(VAR, DATA);

fprintf('\n---------------------------------------------\n');
fprintf('  Newton–Raphson on manifold (LM + Armijo)\n');
fprintf('  Max iterations: %d\n', DATA.NEWTON_RAPHSON.NMAXiter);
fprintf('---------------------------------------------\n');

% -------------------------------------------------------------------------
% 1) Newton iterations
% -------------------------------------------------------------------------
while kiter <= DATA.NEWTON_RAPHSON.NMAXiter
    
    DATA.kiter = kiter;   % expose iter to called kernels if needed
    
    % (Optional) ECM clustered weights retrieval
    if notempty(OPERFE,'wSTs_cluster')
        if ~isfield(OPERFE.wSTs_cluster,'DATA_regress') || ...
                ~isfield(OPERFE.wSTs_cluster.DATA_regress,'IndexLinear_damageMODEL')
            [OPERFE, q_LAT] = RetrieveWeightsMAWECM_plast_LARGE(VAR,OPERFE);
        else
            [OPERFE, q_LAT] = RetrieveWeightsMAWECM_DAMAGE(VAR,OPERFE);
        end
    else
        q_LAT = [];
    end
    
    % ---------------------------------------------------------------------
    % 1.1) Map q -> u via τ(q), assemble residual and project to reduced
    % ---------------------------------------------------------------------
    [VAR, tauNON_q, tauNONder_q, tauNONder2_q, tauDER_all] = ...
        update_state_and_project(VAR, DATA, OPERFE, MATPRO, VARint_n, DOFl);
    
    % External load on reduced space (possibly +hydro below)
    [VAR_FINT_DOFL, VAR_RESID_DOFL, VAR_FEXT_DOFl] = reduced_vectors(VAR, OPERFE, DOFl);
    
    % Early convergence check (cheap)
    [res_ok, step_ok] = robust_stopping(VAR_RESID_DOFL, VAR.DISP(DOFl), norm_step, DATA, VAR_FEXT_DOFl);
    if res_ok && step_ok
        CONVERGED = 1;
    end
    
    % Optional dynamic/hydro contributions (kept as your policy)
    if DATA.ISDYNAMIC == 1 && ~isempty(VAR.RESID)
        error('Dynamic path not adapted yet (intentional).');
    end
    if notempty(OPERFE,'HYDRO')
        error('Hydro path not adapted yet (intentional).');
    end
    
    % Snapshot history
    DATA = UpdateIterativeSnapshot(DATA,VAR,1+kiter);
    
    % If converged, finalize stresses/energies and exit
    if CONVERGED == 1
        VAR = finalize_stresses_and_energies(VAR, DATA);
        break
    end
    
    % ---------------------------------------------------------------------
    % 1.2) Build full tangent and project to reduced; add optional terms
    % ---------------------------------------------------------------------
    Kextend = ...
        GetStiffnessMatrix(OPERFE, DATA, VAR, VAR.FGRADST, VAR.CELASTST); %#ok<ASGLU>
    % (Note: we rely on your kernels to have populated VAR.FGRADST and VAR.CELASTST
    %        inside ResidualFromDisplacementsVARw; if not, pass via outputs as in your code.)
    
    K = tauDER_all' * Kextend * tauDER_all;
    
    % Optional manifold-geometry contribution
    if USE_GEOMETRIC_term_K == 1 && ~isempty(tauNONder2_q)
        K = add_geometric_term(K, tauNONder2_q, VAR.RESID_extend, OPERFE, DOFl);
    end
    
    % Optional ECM weights contribution
    if USE_WEIGHT_TERM_k == 1 && notempty(VAR,'Kred_w_LL')
        Kred_w_LL = VAR.Kred_w_LL;     % passed out by your ResidualFromDisplacementsVARw
        K(DOFl,DOFl) = K(DOFl,DOFl) + Kred_w_LL;
    end
    
    % Optional dynamic/hydro tangent changes (kept as in your codebase)
    if DATA.ISDYNAMIC == 1
        K = KstifflDynamicPart(VAR,OPERFE,DATA.TIME_INTEGRATION_PARAMETERS,DATA.DeltaT,K);
    end
    if notempty(OPERFE,'HYDRO')
        error('Hydro tangent contribution not adapted yet (intentional).');
    end
    
    % ---------------------------------------------------------------------
    % 1.3) LM-regularized step on DOFl and Armijo backtracking
    % ---------------------------------------------------------------------
    r_do_fl = VAR_RESID_DOFL;
    Kred    = K(DOFl,DOFl);
    
   
    
    % LM-regularization
    lam     = DATA.LM_lambda;
    Kreg    = Kred + lam*speye(length(DOFl));
    deltaq  = - solve_safely(Kreg, r_do_fl);
    
    % Merit function at current state
    phi0    = 0.5*(r_do_fl.'*r_do_fl);
    
    % Backtracking with trial re-assembly using your kernels (consistent)
    [accepted, VAR_try, phi_try, alpha, lam_new, dq_eff] = ...
        line_search_armijo(VAR, deltaq, DOFl, phi0, DATA, OPERFE, MATPRO, ...
        VARint_n, tauNONder_q);
    
    % Accept / update damping
    if accepted
        VAR  = VAR_try;
        norm_step     = norm(dq_eff);
        DATA.LM_lambda = max(lam_new, DATA.LM_lambda_min);
    else
        % Try a stronger damping one-shot (no backtracking)
        lam_boost = min(max(10*DATA.LM_lambda, 1e-8), DATA.LM_lambda_max);
        Kreg2     = Kred + lam_boost*speye(length(DOFl));
        deltaq2   = - solve_safely(Kreg2, r_do_fl);
        [accepted2, VAR_try2, phi_try2, alpha2, lam_new2, dq_eff2] = ...
            line_search_armijo(VAR, deltaq2, DOFl, phi0, DATA, OPERFE, MATPRO, ...
            VARint_n, tauNONder_q);
        if accepted2
            VAR  = VAR_try2;
            norm_step      = norm(dq_eff2);
            DATA.LM_lambda = max(lam_new2, DATA.LM_lambda_min);
        else
            % Commit the best seen (no improvement): keep current VAR, increase lambda
            DATA.LM_lambda = min(lam_boost*10, DATA.LM_lambda_max);
            norm_step      = norm(deltaq2);
        end
    end
    
     [VAR, DATA, norm_step] = handle_stagnation( ...
    accepted, DATA, K, DOFl, r_do_fl, VAR, phi0, ...
    OPERFE, MATPRO, VARint_n, tauNONder_q, norm_step);
    
    % ---------------------------------------------------------------------
    % 1.4) Convergence check (residual + step)
    % ---------------------------------------------------------------------
    [VAR_FINT_DOFL, VAR_RESID_DOFL, VAR_FEXT_DOFl] = reduced_vectors(VAR, OPERFE, DOFl);
    [res_ok, step_ok] = robust_stopping(VAR_RESID_DOFL, VAR.DISP(DOFl), norm_step, DATA, VAR_FEXT_DOFl);
    if res_ok && step_ok
        CONVERGED = 1;
        VAR = finalize_stresses_and_energies(VAR, DATA);
        break
    end
    % Diagnostics
    resnorm = norm(VAR.RESID(DOFl));
    fprintf('Iter %2d | ‖r‖ = %.3e | ‖Δq‖ = %.3e | λ = %.2e\n', ...
        kiter, resnorm, norm_step, DATA.LM_lambda);
    
    % Iterate
    kiter = kiter + 1;
end

% If exited without convergence, still report last energies if needed
if CONVERGED ~= 1
    % nothing else to do; caller decides next action
end

end
% =========================================================================
% ===========================  AUXILIARIES  ===============================
% =========================================================================

function s = notempty(S,field)
s = isstruct(S) && isfield(S,field) && ~isempty(S.(field));
end

function DATA = set_default(DATA, name, val)
if ~isfield(DATA,name) || isempty(DATA.(name))
    DATA.(name) = val;
end
end

function VARint_n = pull_internal_variables(VAR, DATA)
VARint_n = struct();
if notempty(DATA,'ListFieldInternalVariables')
    for iINT = 1:numel(DATA.ListFieldInternalVariables)
        nm = DATA.ListFieldInternalVariables{iINT};
        if isfield(VAR,nm), VARint_n.(nm) = VAR.(nm); end
    end
end
end

function [VAR, tauNON_q, tauNONder_q, tauNONder2_q, tauDER_all] = ...
    update_state_and_project(VAR, DATA, OPERFE, MATPRO, VARint_n, DOFl)
% Evaluate τ, build extended displacement, assemble residual, and project.
qL  = VAR.DISP(DOFl);
qR  = VAR.DISP(OPERFE.DOFr);

[tauNON_q, tauNONder_q, tauNONder2_q] = feval( ...
    DATA.DATA_evaluateTAU_and_DER.nameFunctionEvaluate, ...
    qL, DATA.DATA_evaluateTAU_and_DER);

% Build extended displacement u = [τ(qL); qR]
VAR.DISP_q   = VAR.DISP;
VAR.DISP     = [tauNON_q; qR];
VAR.FEXT     = VAR.FEXT_extended;

% Build block τ′ including identity on DOFr
tauDER_all = build_tauDER_all(tauNONder_q, length(OPERFE.DOFr));

% Assemble residuals via your kernel (STATIC)
[VAR, VAR.CELASTST, VAR.FGRADST, detFgrad, Kred_w_LL] = ...
    ResidualFromDisplacementsVARw(OPERFE, VAR, MATPRO, DATA, VARint_n, tauDER_all, DOFl); %#ok<ASGLU>
VAR.Kred_w_LL     = Kred_w_LL;    % store if provided by your kernel
VAR.RESID_extend  = VAR.RESID;

% Project forces/residual to reduced space
VAR.FINT  = tauDER_all' * VAR.FINT;
VAR.FEXT  = tauDER_all' * VAR.FEXT;
VAR.RESID = tauDER_all' * VAR.RESID;

% Restore reduced coordinates vector (qL,qR)
VAR.DISP_EXTENDED = VAR.DISP;
VAR.DISP          = VAR.DISP_q;
end

function tauDER_all = build_tauDER_all(tauNONder_q, nDOFr)
% Build block derivative matrix:
%   [ τ′(qL)    0
%     0         I ]
RR = speye(nDOFr,nDOFr);
LR = sparse(size(tauNONder_q,1), nDOFr);
RL = sparse(nDOFr, size(tauNONder_q,2));
tauDER_all = [tauNONder_q, LR; RL, RR];
end

function [VAR_FINT_DOFL, VAR_RESID_DOFL, VAR_FEXT_DOFl] = reduced_vectors(VAR, OPERFE, DOFl)
if isempty(OPERFE.DOFm)
    VAR_FINT_DOFL  = VAR.FINT(DOFl);
    VAR_RESID_DOFL = VAR.RESID(DOFl);
    VAR_FEXT_DOFl  = VAR.FEXT(DOFl);
else
    VAR_FINT_DOFL  = OPERFE.A' * VAR.FINT;
    VAR_RESID_DOFL = OPERFE.A' * VAR.RESID;
    VAR_FEXT_DOFl  = OPERFE.A' * VAR.FEXT;
end
end

function K = add_geometric_term(K, tauNONder2_q, RESID_ext, OPERFE, DOFl)
% K_geo = sum_i ( τ″_i(q) * r_i(u) ) added to DOFl block
K_geo = 0;
n_general = length(RESID_ext) - length(OPERFE.DOFr);
for i = 1:n_general
    % tauNONder2_q(i,:,:) * RESID_ext(i)  -> [nL x nL] block
    K_geo = K_geo + squeeze(tauNONder2_q(i,:,:)) * RESID_ext(i);
end
K(DOFl,DOFl) = K(DOFl,DOFl) + K_geo;
end

function x = solve_safely(A,b)
% Small helper to safely solve A x = b (prefers backslash; fallback PCG).
try
    x = A \ b;
catch
    % crude fallback
    x = pcg(A, b, 1e-12, 200);
end
end

function [accepted, VAR_out, phi_try, alpha, lam_new, dq_eff] = ...
    line_search_armijo(VAR_in, deltaq, DOFl, phi0, DATA, OPERFE, MATPRO, VARint_n, tauNONder_q_hint)
% Armijo backtracking on φ = 1/2 ||r_DOFl||^2 with LM damping adaptation.
% Reassembles residual consistently via your kernels at trial points.

alpha   = 1.0;
c_arm   = DATA.armijo_c;
rho     = DATA.armijo_rho;
alpha_min = DATA.armijo_alpha_min;

lam_new = max(DATA.LM_lambda*0.3, DATA.LM_lambda_min);  % optimistic decrease upon success

VAR_bak = VAR_in;
phi_try = inf;
accepted = false;
dq_eff = 0;

while alpha >= alpha_min
    % Trial update on reduced coords
    VAR_trial = VAR_bak;
    VAR_trial.DISP(DOFl) = VAR_trial.DISP(DOFl) + alpha*deltaq;
    
    % Re-evaluate τ, u, residual, and project (consistent)
    qL  = VAR_trial.DISP(DOFl);
    qR  = VAR_trial.DISP(OPERFE.DOFr);
    
    %         [tauNON_q, tauNONder_q, ~] = feval( ...
    %             VAR_bak.DATA_evaluateTAU_and_DER.nameFunctionEvaluate, ...
    %             qL, VAR_bak.DATA_evaluateTAU_and_DER);
    
    [tauNON_q, tauNONder_q, ~] = feval( ...
        DATA.DATA_evaluateTAU_and_DER.nameFunctionEvaluate, ...
        qL, DATA.DATA_evaluateTAU_and_DER);
    
    VAR_trial.DISP = [tauNON_q; qR];
    VAR_trial.FEXT = VAR_trial.FEXT_extended;
    
    tauDER_all = build_tauDER_all(tauNONder_q, length(OPERFE.DOFr));
    
    [VAR_trial,~,~,~,~] = ResidualFromDisplacementsVARw( ...
        OPERFE, VAR_trial, MATPRO, DATA, VARint_n, tauDER_all, DOFl);
    
    VAR_trial.RESID_extend = VAR_trial.RESID;
    VAR_trial.FINT = tauDER_all' * VAR_trial.FINT;
    VAR_trial.FEXT = tauDER_all' * VAR_trial.FEXT;
    VAR_trial.RESID = tauDER_all' * VAR_trial.RESID;
    
    % Reduced residual and merit
    if isempty(OPERFE.DOFm)
        r_try = VAR_trial.RESID(DOFl);
    else
        r_try = OPERFE.A' * VAR_trial.RESID;
    end
    phi_try = 0.5*(r_try.'*r_try);
    
    % Armijo: φ(x+αΔ) <= φ(x) - c α ||r||^2
    if phi_try <= (phi0 - c_arm*alpha*(2*phi0))
        accepted = true;
        VAR_out  = VAR_trial;
        dq_eff   = alpha*deltaq;
        return
    end
    
    % Backtrack
    alpha = rho*alpha;
end

VAR_out  = VAR_in;
lam_new  = min(DATA.LM_lambda*10, DATA.LM_lambda_max); % increase damping if no step accepted
dq_eff   = 0;
end

function [res_ok, step_ok] = robust_stopping(r_DOFl, q_DOFl, norm_step, DATA, fext_DOFl)
% Residual + step scaled tests
if nargin < 5 || isempty(fext_DOFl), fext_DOFl = 0; end
res_tol  = DATA.TOL_RES * (1 + norm(fext_DOFl,2));
step_tol = DATA.TOL_STEP * (1 + norm(q_DOFl,2));
res_ok   = (norm(r_DOFl,2) <= res_tol);
step_ok  = (norm_step <= step_tol);
end

function VAR = finalize_stresses_and_energies(VAR, DATA)
% Postprocess stresses and energies (mirrors your original code paths)
if DATA.SMALL_STRAIN_KINEMATICS == 0
    PK1_forCauchy = VAR.PK1STRESS;
else
    PK1_forCauchy = VAR.PK2STRESS;
end
[VAR.VONMISES_CAUCHY_STRESS, VAR.CAUCHY_STRESS] = ...
    VonMisesCauchyStresses(PK1_forCauchy, VAR.FGRADST, DATA.MESH.ndim, VAR.detFgrad, DATA);

if DATA.ISDYNAMIC == 1
    VAR = UpdateVelocityAcceleration(VAR,DATA.TIME_INTEGRATION_PARAMETERS,DATA.DeltaT);
end
VAR = KineticAndStrainEnergyGet_mnfold(VAR,OPERFE_placeholder(),DATA,VAR.FGRADST); %#ok<NASGU>
end

function O = OPERFE_placeholder()
% Placeholder to satisfy signature of KineticAndStrainEnergyGet_mnfold if needed;
% replace with actual OPERFE if your function requires it.
O = struct();
end

function [VAR, DATA, norm_step] = handle_stagnation( ...
    accepted, DATA, K, DOFl, r_do_fl, VAR, phi0, ...
    OPERFE, MATPRO, VARint_n, tauNONder_q, norm_step)
%HANDLE_STAGNATION  Breaks LM saturation / stagnation.
%
% Called when either the step was rejected or lambda hit its upper bound.
% Tries (A) bold low-lambda step, (B) Broyden-good secant update.
%
% INPUTS
%   accepted    : flag from previous Armijo call
%   DATA, K, DOFl, r_do_fl, VAR, phi0, OPERFE, MATPRO, VARint_n, tauNONder_q
%   norm_step   : previous step norm
%
% OUTPUTS
%   VAR, DATA   : possibly updated after rescue step
%   norm_step   : updated step norm
%
% -------------------------------------------------------------------------

if (~accepted) || (DATA.LM_lambda >= DATA.LM_lambda_max - eps)
    % Count consecutive stagnation
    if ~isfield(DATA,'_stagn_iters'), DATA.stagn_iters = 0; end
    DATA.stagn_iters = DATA.stagn_iters + 1;
    
    M = 5;  % how many consecutive stagnant iters trigger rescue
    if DATA.stagn_iters >= M
        % ===== STRATEGY A: bold step (reset lambda way down) =====
        DATA.LM_lambda = max(1e-10, DATA.LM_lambda_min);
        Kred = K(DOFl,DOFl);
        Kreg = Kred + DATA.LM_lambda*speye(length(DOFl));
        deltaq_bold = - solve_safely(Kreg, r_do_fl);
        
        armijo_c_backup = DATA.armijo_c;
        DATA.armijo_c   = 1e-6;   % relax Armijo for this rescue
        [acc_bold, VAR_bold, phi_bold, alpha_bold, lam_new_bold, dq_eff_bold] = ...
            line_search_armijo(VAR, deltaq_bold, DOFl, phi0, DATA, ...
            OPERFE, MATPRO, VARint_n, tauNONder_q);
        DATA.armijo_c = armijo_c_backup;
        
        if acc_bold
            VAR = VAR_bold;
            norm_step = norm(dq_eff_bold);
            DATA.LM_lambda = max(lam_new_bold, DATA.LM_lambda_min);
            DATA.stagn_iters = 0;
            return;
        end
        
        % ===== STRATEGY B: Broyden-good update on reduced space =====
        if ~isfield(DATA,'_use_broyden') || DATA.use_broyden==0
            DATA.use_broyden = 1;
            DATA.B = Kred;   % initialise secant matrix
        end
        
        B = DATA.B;
        dq_b = - solve_safely(B, r_do_fl);
        
        [acc_b, VAR_b, phi_b, alpha_b, lam_dummy, dq_eff_b] = ...
            line_search_armijo(VAR, dq_b, DOFl, phi0, DATA, ...
            OPERFE, MATPRO, VARint_n, tauNONder_q);
        
        if acc_b
            % Secant update: B_{k+1} = B_k + ((y - B_k s) s^T)/(s^T s)
            r_before = r_do_fl;
            [~, r_after, ~] = reduced_vectors(VAR_b, OPERFE, DOFl);
            s = alpha_b * dq_b;
            y = r_after - r_before;
            denom = s.'*s + eps;
            DATA.B = B + ((y - B*s) * (s.')) / denom;
            
            VAR = VAR_b;
            norm_step = norm(dq_eff_b);
            DATA.stagn_iters = 0;
        else
            DATA.LM_lambda = min(10*max(DATA.LM_lambda,1e-8), DATA.LM_lambda_max);
        end
     
else
    % Progress made → reset counter
    if isfield(DATA,'_stagn_iters'), DATA.stagn_iters = 0; end
end
end

