function [VAR,CONVERGED,DATA] = NewtonRapshonStatic_MHROM_LSearch_GROK_disab(DATA,OPERFE,VAR,MATPRO,DOFl)
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

% ToleranceNR_OLD_REL = DATA.NEWTON_RAPHSON.TOL_FORCES_REL ;
% ToleranceNR_OLD_ABS = DATA.NEWTON_RAPHSON.TOL_DISPLAC_ABS ;
%
%

while kiter <= DATA.NEWTON_RAPHSON.NMAXiter
    
    [FintL,Rl,VAR,OPERFE,q_LAT,tauNON_q,tauNONder_q,tauNONder2_q,celastST,FgradST...
        ,detFgrad,Kred_w_LL,DATA,VAR_FEXT_DOFl,tauNONallDER_q]...
        = GetResidual_MHROM(OPERFE,VAR,VARint_n,DATA,MATPRO,DOFl,kiter) ;
    
    if isempty(celastST) && DATA.INTERNAL_FORCES_USING_precomputed_Kstiff == 0
        disp('Not converged ..., because of Negative JAcobians (you may disable this by setting PROCEED_WITH_NEGATIVE_JACOBIANS=1)')
        DATA = ReshapeIterativeSnapshot(DATA,VAR,kiter) ;
        break
    else
        % 7. Check convergence
        [norm_Rl, ~,CONVERGED]= CheckConvergenceLSTRq(FintL,VAR_FEXT_DOFl,Rl,kiter,normd,DATA,q_LAT) ;
        %
        %         DATA.NEWTON_RAPHSON.TOL_FORCES_REL  =  ToleranceNR_OLD_REL  ;
        %        DATA.NEWTON_RAPHSON.TOL_DISPLAC_ABS =  ToleranceNR_OLD_ABS   ;
        %
        
        
        if CONVERGED == 1
            VAR = UpdateVariablesNR(DATA,VAR,FgradST,detFgrad,OPERFE) ;
            break
        else
            [Kll,KllGEO,KllW] = GetJacobianMatrix_MHROM(OPERFE, DATA, VAR, FgradST, celastST,tauNONallDER_q,tauNONder2_q,DOFl,Kred_w_LL) ;
            
            % Parameters for line search and regularization
            %   max_reg_tries = 10;
            max_ls_tries = 10;
            c_ls = 1e-4;
            beta_ls = 0.5;
            reg_lambda = 0;
            eps_reg =  1e-6;
            NumberOfTrials_lambda = 10 ;
            Max_lambda_regularization = 1e-2;
            converged_step = false;
            
            m0 = 0.5 * (Rl' * Rl);
            norm_Rl_before = norm(Rl) ; 
            e = eig(Kll);
            % Minimum eigenvalue
            
            min_re = min(real(e));
            
            
            %
            %max_re = max(real(e));
%             if min_re > eps_reg
%                 % Scenario  1: No need for regularization (symmetric part is positive definite)
%                 lambdaVECTOR = 0 ;
%             else
%                
%                   e_noGEO = real(eig(Kll_noGEO));
%                   e_noW = real(eig(Kll_noW));                
%                   e_MAT = real(eig(Kll_mat));
% %                 
% %                 lambdaVECTOR =abs(min_re) +  linspace(eps_reg,Max_lambda_regularization,NumberOfTrials_lambda) ;
%             end
%             

                 Kll_noGEO = Kll-KllGEO ; 
                 Kll_noW = Kll-KllW ;
                Kll_mat = Kll-KllW-KllGEO ;
                
                

            MatricesToTry =   {Kll_mat} ; %{Kll,Kll_noGEO,Kll_noW,Kll_mat} ; 




            for iMATRIX = 1:length(MatricesToTry)
                % Compute eigenvalues to check for regularization
                
                
                %    current_lambda = reg_lambda + max(0, -min_re + eps_reg);
              %  current_lambda =   lambdaVECTOR(reg_try) ;
                % Regularized Jacobian
              %  Kreg = Kll + current_lambda * eye(size(Kll, 1));
                Kreg = MatricesToTry{iMATRIX} ; 
                % Search direction
                %  try
                delta_dL = -Kreg \ Rl;
                %                 catch
                %                     % If still singular, increase regularization and continue
                %                     if reg_lambda == 0
                %                         reg_lambda = eps_reg;
                %                     else
                %                         reg_lambda = reg_lambda * 2;
                %                     end
                %                     continue;
                % end
                
                % Line search with backtracking (Armijo condition)
                alpha = 1;
                converged_ls = false;
                for ls_try = 1:max_ls_tries
                    VAR_trial = VAR;
                    VAR_trial.DISP(DOFl) = VAR.DISP(DOFl) + alpha * delta_dL;
                    
                    if ~isempty(OPERFE.DOFm)
                        VAR_trial.DISP(OPERFE.DOFr) = OPERFE.G * VAR_trial.DISP(OPERFE.DOFm) + OPERFE.uBAR;
                    end
                    
                    [~, Rl_trial, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] ...
                        = GetResidual_MHROM(OPERFE, VAR_trial, VARint_n, DATA, MATPRO, DOFl, kiter);
                    
                    m_trial = 0.5 * (Rl_trial' * Rl_trial);
                     
                    
                    if (m_trial <= m0 + alpha * c_ls * (-2 * m0)) ||  norm(Rl_trial) <  norm_Rl_before
                        converged_ls = true;
                        break;
                    end
                    
                    alpha = alpha * beta_ls;
                end
                
                if converged_ls
                    VAR = VAR_trial;
                    normd = norm(alpha * delta_dL);
                    converged_step = true;
                    if alpha ~= 1  ||   iMATRIX ~=1
                        fprintf('LS converged: 1/alpha = %.d, Matrix Kll option = %.d \n', 1/alpha,iMATRIX);
                    end
                    
                    break;
                else
                    %   fprintf('Line Search NOT converged with lambda_reg = %.3e\n',  current_lambda);
                    % Increase progressive regularization for next try
                    %                     if reg_lambda == 0
                    %                         reg_lambda = eps_reg;
                    %                     else
                    %                         reg_lambda = reg_lambda * 2;
                    %
                    %                     end
                    disp('')
                end
            end
            
            if ~converged_step
                disp('Newton-Raphson step failed after max regularization and line search attempts.');
                DATA = ReshapeIterativeSnapshot(DATA, VAR, kiter);
                break
                %                 Tolmin = 1e-4;
                %                 if norm_Rl < Tolmin
                %                     % We pretend it has converged....
                %                     TOLlocNR = Tolmin;
                %
                %                     DATA.NEWTON_RAPHSON.TOL_FORCES_REL  = TOLlocNR ;
                %        DATA.NEWTON_RAPHSON.TOL_DISPLAC_ABS =  TOLlocNR ;
                %
                % %               ToleranceNR_OLD_REL = DATA.NEWTON_RAPHSON.TOL_FORCES_REL ;
                % % ToleranceNR_OLD_ABS = DATA.NEWTON_RAPHSON.TOL_DISPLAC_ABS ;
                % %
                %
                %
                %                 else
                %                     break
                %                 end
                
            end
            
        end
        kiter = kiter + 1;
        
    end
end