function [FintL,Rl,VAR,OPERFE,q_LAT,tauNON_q,tauNONder_q,tauNONder2_q,celastST,FgradST,...
    detFgrad,Kred_w_LL,DATA,VAR_FEXT_DOFl,tauNONallDER_q]...
    = GetResidual_MHROM(OPERFE,VAR,VARint_n,DATA,MATPRO,DOFl,kiter)
%--------------------------------------------------------------------------
% GetResidual_MHROM
% -----------------
% PURPOSE
%   Assemble the reduced residual, internal forces and auxiliary tensors for a
%   manifold-based HROM at a given Newton iteration. The routine:
%     1) (Optional) updates ECM weights when MAW-ECM clustering is enabled,
%     2) evaluates the manifold map τ(qL) and its derivatives τ′(qL), τ″(qL),
%     3) builds the extended reduced state [τ(qL); qR],
%     4) calls the FE/HROM residual kernel at the current iterate,
%     5) projects forces/residuals to the tangent manifold via τ′(qL),
%     6) returns DOFl blocks and bookkeeping needed by the Newton solve.
%
% CONTEXT
%   • Manifold ROM with tangent-space projection u = [τ(qL); qR].
%   • Hyperreduction via ECM, with optional **latent-dependent** weight
%     scheduling (MAW-ECM) handled through OPERFE.wSTs_cluster.
%
% INPUTS
%   OPERFE : Operators/partitions and ECM data
%            - DOFr                 : indices of constrained reduced DOFs
%            - wSTs_cluster (opt.)  : struct enabling MAW-ECM scheduling
%               · DATA_regress.IndexLinear_damageMODEL (flag for damage path)
%            - (other fields used by ResidualFromDisplacementsVARw)
%   VAR    : State at entry
%            - DISP           : current reduced coordinates [qL;qR]
%            - FEXT_extended  : external forces in extended space
%            - (internal variables possibly carried in VARint_n)
%   VARint_n : Internal variables at previous (time) step (structure of arrays)
%   DATA   : Simulation controls and τ-evaluator descriptor
%            - DATA_evaluateTAU_and_DER.nameFunctionEvaluate : function handle
%            - DATA_evaluateTAU_and_DER : packed data for τ, τ′, τ″ evaluation
%            - ISDYNAMIC            : if 1 → dynamic residual (not adapted here)
%            - kiter                : current Newton iteration counter (set here)
%   MATPRO : Material parameters used by the residual kernel
%   DOFl   : Indices of free reduced DOFs (qL block)
%   kiter  : Current Newton iteration (integer)
%
% OUTPUTS
%   FintL             : Reduced internal force vector on DOFl
%   Rl                : Reduced residual on DOFl
%   VAR               : Updated state (includes projected FINT/FEXT/RESID and caches)
%   OPERFE            : Possibly updated (if MAW-ECM reweights were retrieved)
%   q_LAT             : Latent abscissa used for MAW-ECM selection (if active)
%   tauNON_q          : τ(qL)          (displacement embedding)
%   tauNONder_q       : τ′(qL)         (Jacobian of embedding)
%   tauNONder2_q      : τ″(qL)         (second derivative, if provided)
%   celastST, FgradST : Constitutive/kinematic tensors from residual kernel
%   detFgrad          : Determinant(s) of deformation gradient (for checks)
%   Kred_w_LL         : ECM-weighted reduced tangent block (LL)
%   DATA              : Updated with iterative snapshot via UpdateIterativeSnapshot
%   VAR_FEXT_DOFl     : Reduced external force on DOFl (incl. pressure if used)
%   tauNONallDER_q    : Block-Jacobian for extended coordinates [τ′  0; 0  I]
%
% EXECUTION FLOW
%   (A) MAW-ECM weight scheduling (optional)
%       if OPERFE.wSTs_cluster:
%           - RetrieveWeightsMAWECM_plast_LARGE or ..._DAMAGE → update OPERFE, q_LAT.
%   (B) Manifold evaluation
%       - qL = VAR.DISP(DOFl); qR = VAR.DISP(OPERFE.DOFr)
%       - [τ, τ′, τ″] = feval(DATA.DATA_evaluateTAU_and_DER.nameFunctionEvaluate, qL, ...)
%       - Build extended reduced state: VAR.DISP = [τ; qR]
%       - Pack τ′ into block-Jacobian τall′ = [τ′ 0; 0 I_r]
%   (C) Residual kernel
%       - DATA.kiter = kiter; VAR.FEXT = VAR.FEXT_extended
%       - [VAR, celastST, FgradST, detFgrad, Kred_w_LL] =
%             ResidualFromDisplacementsVARw(OPERFE, VAR, MATPRO, DATA, VARint_n, τall′, DOFl)
%   (D) Tangent-space projection
%       - VAR.FINT = τall′ᵀ VAR.FINT;  VAR.FEXT = τall′ᵀ VAR.FEXT;  VAR.RESID = τall′ᵀ VAR.RESID
%       - Restore VAR.DISP back to reduced coordinates (undo extended state)
%   (E) Options not adapted here
%       - Dynamic residual path (ISDYNAMIC == 1) → error
%       - Hydrostatic follower loads (OPERFE.HYDRO) → error
%   (F) External loads on DOFl
%       - VAR_FEXT_DOFl = VAR.FEXT(DOFl) (or + pressure term if enabled)
%   (G) DOFl block extraction
%       - If no DOFm: FintL = VAR.FINT(DOFl); Rl = VAR.RESID(DOFl)
%       - Else (affine/periodic BCs): left-multiply by Aᵀ
%   (H) Snapshot
%       - DATA = UpdateIterativeSnapshot(DATA, VAR, 1 + kiter)
%
% NOTES / SAFEGUARDS
%   • If detFgrad is empty and INTERNAL_FORCES_USING_precomputed_Kstiff == 0,
%     upper layers typically flag “negative Jacobians” and stop the solve.
%   • τ″(qL) may be empty; downstream geometric terms should check before using it.
%   • This function does not perform line search or regularization; those are
%     responsibilities of the Newton driver.
%
% DEPENDENCIES
%   RetrieveWeightsMAWECM_plast_LARGE, RetrieveWeightsMAWECM_DAMAGE,
%   ResidualFromDisplacementsVARw, UpdateIterativeSnapshot.
%
% AUTHOR
%   Joaquín A. Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%
% COMMENTS UPDATED BY: ChatGPT-5
%   Date of update: 07-NOV-2025, Barcelona, Spain
%--------------------------------------------------------------------------

if  ~isempty(OPERFE.wSTs_cluster)
    if ~isfield(OPERFE.wSTs_cluster.DATA_regress,'IndexLinear_damageMODEL')
        [OPERFE,q_LAT] = RetrieveWeightsMAWECM_plast_LARGE(VAR,OPERFE) ;
    else
        [OPERFE,q_LAT] = RetrieveWeightsMAWECM_DAMAGE(VAR,OPERFE) ;
    end
end


% Residual forces from displacements  (STATIC )
%qALL = VAR.DISP ;
qL = VAR.DISP(DOFl) ;
qR = VAR.DISP(OPERFE.DOFr) ;
VAR.DISP_q = VAR.DISP ; % these are the actual reduced coordinates
% Now we evaluate the coefficientes of the modal expansion
[tauNON_q,tauNONder_q,tauNONder2_q] = feval(DATA.DATA_evaluateTAU_and_DER.nameFunctionEvaluate,qL,DATA.DATA_evaluateTAU_and_DER);

% Now we create the displacement using "extended" coordinate
VAR.DISP = [tauNON_q;qR] ;
% Now we can invoke the function that computes the residual
DATA.kiter = kiter;
VAR.FEXT = VAR.FEXT_extended;

tauNONallDER_q_RR = speye(length(OPERFE.DOFr),length(OPERFE.DOFr)) ;
tauNONallDER_q_LR = sparse(size(tauNONder_q,1),length(OPERFE.DOFr));
tauNONallDER_q_RL =  sparse(length(OPERFE.DOFr),size(tauNONder_q,2));

tauNONallDER_q = [tauNONder_q,tauNONallDER_q_LR
    tauNONallDER_q_RL, tauNONallDER_q_RR ] ;
[VAR,celastST,FgradST,detFgrad,Kred_w_LL] =  ResidualFromDisplacementsVARw(OPERFE,VAR,MATPRO,DATA,VARint_n,tauNONallDER_q,DOFl) ;
% Now we have to project some of these variables on the tangent manifold
% To this end, we need tauNONallDER_q
% This is a matrix of size nmodes_extended x ndof


% Projection onto the reduced coefficients
VAR.RESID_extend = VAR.RESID ;
VAR.FINT = tauNONallDER_q'*VAR.FINT ;
VAR.FEXT = tauNONallDER_q'*VAR.FEXT ;
VAR.RESID = tauNONallDER_q'*VAR.RESID ;
% Go back to the reduced coordinates
VAR.DISP_EXTENDED = VAR.DISP ;
VAR.DISP = VAR.DISP_q ;

if DATA.ISDYNAMIC == 1 && ~isempty(VAR.RESID)
    % Dynamic part of the residual
    error('Option not adapted yet, 2-July-2025')
    VAR =  ResidualDynamicPart(VAR,OPERFE,DATA.TIME_INTEGRATION_PARAMETERS,DATA.DeltaT) ;
end
if ~isempty(OPERFE.HYDRO)
    error('Option not adapted yet, 2-July-2025')
    % Hydrostatic forces
    % \FextPR{}{}   = - \NbST{T}  \wSTft{}  \tPRESSst}
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/DynamicFloatingBodies.tex
    % ------------------------
    [VAR.FEXTpress,VAR.tPRESS,kPRESSst ]= PressureForcesFromDisplacements(DATA,OPERFE,VAR) ;
    VAR.RESID = VAR.RESID - VAR.FEXTpress ;
    VAR_FEXT_DOFl = VAR.FEXT(DOFl) + VAR.FEXTpress(DOFl);
else
    VAR_FEXT_DOFl = VAR.FEXT(DOFl) ;
end

% Update snapshots of iterative results
DATA = UpdateIterativeSnapshot(DATA,VAR,1+kiter) ;

if isempty(OPERFE.DOFm)
    FintL = VAR.FINT(DOFl) ;
    Rl = VAR.RESID(DOFl) ;
else
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/10_PeriodicQUADL.mlx
    FintL = OPERFE.A'*VAR.FINT ;
    Rl = OPERFE.A'*VAR.RESID ;
end
