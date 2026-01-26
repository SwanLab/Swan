function [VAR,CONVERGED,DATA] = NewtonRapshonStatic_MHROM_LSearch_GROK(DATA,OPERFE,VAR,MATPRO,DOFl)
%--------------------------------------------------------------------------
% NewtonRapshonStatic_MHROM_LSearch_GROK
% --------------------------------------
% PURPOSE
%   Static Newton–Raphson solver on a nonlinear τ(q)-manifold with
%   backtracking line search and spectral regularization of the reduced
%   Jacobian. Designed for Manifold-HROM (multi-latent q) with optional
%   ECM weight scheduling handled elsewhere.
%
% ALGORITHM OVERVIEW
%   Per Newton iteration k:
%     1) Residual & tangent ingredients
%        - GetResidual_MHROM(...) → Fint_L, R_L (on DOFl), τ, τ′, τ″ (if any),
%          celastST, FgradST, and weighted reduced tangent block Kred_w_LL.
%     2) Convergence check
%        - CheckConvergenceLSTRq(...) on DOFl; exit on success.
%     3) Reduced Jacobian
%        - GetJacobianMatrix_MHROM(...) builds Kll (+ optional geometric KllGEO
%          if τ″·r enabled).
%     4) Spectral regularization
%        - If min real eig(Kll) ≤ ε, try a ladder of λ: Kreg = Kll + λ I.
%          Solve δqL = −Kreg \ R_L for each λ candidate.
%     5) Backtracking line search (Armijo)
%        - Merit m(α)=½‖R_L(qL+αδqL)‖²; choose α by α←α·β until
%          m(α) ≤ m(0) + α·c·(−2 m(0)). Accept first (λ,α) that satisfies it.
%     6) Accept & update
%        - Update qL, record step norm; on final convergence, call UpdateVariablesNR.
%   If no (λ,α) pair succeeds, the step fails; snapshot is stored for diagnostics.
%
% INPUTS (summary)
%   DATA   : Solver settings (NEWTON_RAPHSON.*, tolerances, flags such as
%            PROCEED_WITH_NEGATIVE_JACOBIANS, geometric/weight terms, etc.).
%   OPERFE : Operators/partitions/maps (DOFl/DOFr/DOFm, A, G, uBAR), optional hydro.
%   VAR    : State container (DISP, FEXT_extended, internal variables).
%   MATPRO : Constitutive parameters.
%   DOFl   : Indices of free latent coordinates.
%
% OUTPUTS
%   VAR       : Updated state at exit (stresses/energies consistent on success).
%   CONVERGED : 1 if converged within NMAXiter; 0 otherwise.
%   DATA      : Iteration history/snapshots updated for analysis.
%
% KEY FEATURES
%   • Multi-latent qL on τ(q) manifold.
%   • Armijo backtracking line search in reduced space.
%   • Spectral (Tikhonov-style) regularization Kll + λ I with adaptive λ grid.
%   • Optional geometric stiffness via τ″·r in the reduced tangent.
%   • Compatible with ECM-weighted tangents returned by residual assembly.
%
% CONVERGENCE & SAFETY
%   • Uses CheckConvergenceLSTRq on DOFl.
%   • If negative Jacobians detected and disallowed, exits non-converged
%     after snapshotting the iterate.
%   • Dynamic terms intentionally disabled (static solver).
%
% PRACTICAL TIPS
%   • Keep GetResidual_MHROM and GetJacobianMatrix_MHROM consistent (same weights
%     and linearization point) to maintain descent.
%   • Frequent large λ suggests trying geometric term, refining τ offline, or
%     revisiting ECM weights.
%   • Tune line-search constants (c, β) for robustness vs. speed.
%
% DEPENDENCIES
%   GetResidual_MHROM, GetJacobianMatrix_MHROM, CheckConvergenceLSTRq,
%   UpdateVariablesNR, ReshapeIterativeSnapshot.
%
% HISTORY
%   • 02-Jul-2025  Original single-latent version (Barcelona).
%   • 11-Aug-2025  FST variant with fast τ evaluator (Cartagena).
%   • 13-Aug-2025  MD version (multi-latent q) and ECM auto-switch.
%   • 29-Aug-2025  Comment overhaul (ChatGPT).
%   • 29-Oct-2025  Line-search variant (GROK) with spectral regularization
%                  (UPC Terrassa).
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
            [Kll,KllGEO] = GetJacobianMatrix_MHROM(OPERFE, DATA, VAR, FgradST, celastST,tauNONallDER_q,tauNONder2_q,DOFl,Kred_w_LL) ;
            %Kll = (Kll-KllW) + 0.5*(KllW+KllW');  
            
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
            e = eig(Kll);
            % Minimum eigenvalue 
              
            min_re = min(real(e));
            % Now we contemplate three scenarios
            % 
            %max_re = max(real(e));
            if min_re > eps_reg 
                % Scenario  1: No need for regularization (symmetric part is positive definite)
                lambdaVECTOR = 0 ; 
            else
                lambdaVECTOR =abs(min_re) +  linspace(eps_reg,Max_lambda_regularization,NumberOfTrials_lambda) ; 
            end
            
            for reg_try = 1:length(lambdaVECTOR)
                % Compute eigenvalues to check for regularization
                
                
            %    current_lambda = reg_lambda + max(0, -min_re + eps_reg);
               current_lambda =   lambdaVECTOR(reg_try) ;  
                % Regularized Jacobian
                Kreg = Kll + current_lambda * eye(size(Kll, 1));
                
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
                    
                    if m_trial <= m0 + alpha * c_ls * (-2 * m0)
                        converged_ls = true;
                        break;
                    end
                    
                    alpha = alpha * beta_ls;
                end
                
                if converged_ls
                    VAR = VAR_trial;
                    normd = norm(alpha * delta_dL);
                    converged_step = true;
                    if alpha ~= 1 || reg_lambda ~=0 
                    fprintf('LS converged: 1/alpha = %.d, eigMIN = %.3e,  lambda_reg = %.3e\n', 1/alpha, min_re, current_lambda);
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