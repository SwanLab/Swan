function [VAR,CONVERGED,DATA] = NewtonRapshonStaticLarge_manifoldMDw_fixed(DATA,OPERFE,VAR,MATPRO,DOFl,USE_GEOMETRIC_term_K,...
    USE_WEIGHT_TERM_k)
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

CHECK_convergence =0;

if nargin == 0
    load('tmp1.mat')
    CHECK_convergence =0;
    USE_GEOMETRIC_term_K = 0 ;
    USE_WEIGHT_TERM_k = 0 ;
    %USE_WEIGHT_TERM_k = 0 ;
end



kiter = 1;
CONVERGED = 0 ;
normd = 1e20 ;
% Internal variables (value at previous time step)
% ----------------------
VARint_n = [];
if ~isempty(DATA.ListFieldInternalVariables)
    for iINTVAR = 1:length(DATA.ListFieldInternalVariables)
        NameIntVar= DATA.ListFieldInternalVariables{iINTVAR} ;
        VARint_n.(NameIntVar) = VAR.(NameIntVar) ;
    end
end


% if  ~isempty(OPERFE.wSTs_cluster)
%     q = VAR.DISP(OPERFE.wSTs_cluster.IndexDOFl_q);
%
%
%     [~, idx] = min(abs(OPERFE.wSTs_cluster.q - q));
%     OPERFE.wSTs = OPERFE.wSTs_cluster.Values(:,idx) ;
%     disp('---------------')
%     disp(['Local indexes points nonzero weights =',num2str(find(OPERFE.wSTs~=0)')])
%     disp('---------------')
% end

iterCHANGE_ECM_rule= 1;
nmax_iterCHANGE_ECM_rule =10 ;



while  kiter <=DATA.NEWTON_RAPHSON.NMAXiter
    
    
    if  ~isempty(OPERFE.wSTs_cluster)
        if ~isfield(OPERFE.wSTs_cluster.DATA_regress,'IndexLinear_damageMODEL')
       [OPERFE,q_LAT] = RetrieveWeightsMAWECM_plast_LARGE(VAR,OPERFE) ; 
        else
             [OPERFE,q_LAT] = RetrieveWeightsMAWECM_DAMAGE(VAR,OPERFE) ; 
        end
    end
    
    
    % Residual forces from displacements  (STATIC )
    % Before 9-Feb-2022
    %     [VAR.GLSTRAINS,VAR.PK2STRESS,celastST,VAR.RESID,Fint,FgradST,VAR.PK1STRESS,detFgrad] =...
    %         ResidualFromDisplacements(OPERFE,VAR.DISP,MATPRO,DATA,VAR.FEXT) ;%
    % After 9-Feb-2022
    
    % EVALUATION NONLINEAR FUNCTION RELATING AMPLITUDED DISPLACEMENT MODES
    % WITH THE REDUCED COORDINATES
    %qALL = VAR.DISP ;
    qL = VAR.DISP(DOFl) ;
    qR = VAR.DISP(OPERFE.DOFr) ;
    VAR.DISP_q = VAR.DISP ; % these are the actual reduced coordinates
    % Now we evaluate the coefficientes of the modal expansion
    [tauNON_q,tauNONder_q,tauNONder2_q] = feval(DATA.DATA_evaluateTAU_and_DER.nameFunctionEvaluate,qL,DATA.DATA_evaluateTAU_and_DER);
    
    %     tauNON_q= OPERFE.tauNON(qL);
    %     tauNONder_q = OPERFE.tauNONder(qL);
    %     if ~isempty( OPERFE.tauNONder2)
    %         tauNONder2_q = OPERFE.tauNONder2(qL);
    %         % this is only for the unconstrained DOFs
    %     else
    %         tauNONder2_q = [] ;
    %     end
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
    
    if isempty(celastST) && DATA.INTERNAL_FORCES_USING_precomputed_Kstiff == 0
        disp('Not converged ..., because of Negative JAcobians (you may disable this by setting PROCEED_WITH_NEGATIVE_JACOBIANS=1)')
        DATA = ReshapeIterativeSnapshot(DATA,VAR,kiter) ;
        break
    else
        % 7. Check convergence
        if isempty(OPERFE.DOFm)
            VAR_FINT_DOFL = VAR.FINT(DOFl) ;
            VAR_RESID_DOFL = VAR.RESID(DOFl) ;
        else
            % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/10_PeriodicQUADL.mlx
            VAR_FINT_DOFL = OPERFE.A'*VAR.FINT ;
            VAR_RESID_DOFL = OPERFE.A'*VAR.RESID ;
        end
        [~, ~,CONVERGED]=  CheckConvergenceLSTRq(VAR_FINT_DOFL,VAR_FEXT_DOFl,VAR_RESID_DOFL,kiter,normd,DATA,q_LAT) ;
        
        
        if CONVERGED == 1
            % Compute Von Mises Stress
            if DATA.SMALL_STRAIN_KINEMATICS == 0
                PK1STRESS_forCauchy  = VAR.PK1STRESS ;
            else
                PK1STRESS_forCauchy  = VAR.PK2STRESS ;
            end
            [ VAR.VONMISES_CAUCHY_STRESS,VAR.CAUCHY_STRESS ]= ...
                VonMisesCauchyStresses(PK1STRESS_forCauchy,FgradST,DATA.MESH.ndim,detFgrad,DATA) ;
            if DATA.ISDYNAMIC == 1
                
                % Update velocities and accelerations
                VAR =  UpdateVelocityAcceleration(VAR,DATA.TIME_INTEGRATION_PARAMETERS,DATA.DeltaT) ;
                %  VAR = KineticAndStrainEnergyGet(VAR,OPERFE,DATA,MATPRO) ;
                %   6-May-2023
                %
            end
            VAR = KineticAndStrainEnergyGet_mnfold(VAR,OPERFE,DATA,FgradST) ;  % 6-May-2023
            break
        else
            
            %--------------------------------------------------------------------------
            % 1. Compute full-space tangent stiffness matrix based on current deformation state
            %    This includes material and geometric contributions in the full displacement space.
            Kextend = GetStiffnessMatrix(OPERFE, DATA, VAR, FgradST, celastST) ;
            
            %--------------------------------------------------------------------------
            % 2. Project the full stiffness matrix onto the tangent space of the nonlinear manifold
            %    Following Eq. (184):   K_ROM = τ′(q)^T · K · τ′(q)
            %    This yields the reduced tangent matrix governing the update of q in the ROM.
            
            if kiter== 1
            K = tauNONallDER_q' * Kextend * tauNONallDER_q ;
            
            else
                
            end
            
            
            if CHECK_convergence == 1
                warning('Borrar esto ! ')
                eigBEFORE_2nd = eig(K(DOFl,DOFl))
                
                if any(eigBEFORE_2nd<0)
                    disp('Negative eigenvalue')
                end
            end
            
            %--------------------------------------------------------------------------
            % 3. Optional: Add second-order geometric correction to the stiffness matrix
            %    If τ″(q) is available, use Eq. (185) from the reference to include
            %    geometric stiffness arising from the curvature of the manifold.
            %    This term is: ∑_i [ τ″_i(q) * r_i(u) ]  (scalar-weighted Hessian contributions)
            %                tauNONder2_q = [] ;
            %                disp('borrar esta parte, en la que pongo vacio tauNONder2_q')
            %   tauNONder2_q = [] ;
            if USE_GEOMETRIC_term_K == 1
                K_geo =0 ;
                
                % Loop over generalized DOFs (excluding prescribed DOFs)
                for imodesEXT = 1:(length(VAR.RESID_extend) - length(OPERFE.DOFr))
                    % Geometric stiffness contribution: second derivative * residual component
                    K_geo = K_geo + squeeze(tauNONder2_q(imodesEXT,:,:) * VAR.RESID_extend(imodesEXT));
                end
                
                % Add geometric stiffness contribution to the reduced tangent matrix
                % Only affecting the reduced DOFs (DOFl)
                K(DOFl, DOFl) = K(DOFl, DOFl) + K_geo ;
                
            end
            
            if CHECK_convergence == 1
                warning('Borrar esto ! ')
                eigAFTER = eig(K(DOFl,DOFl))
                
                if any(eigAFTER<0)
                    disp('Negative eigenvalue after geometric contribution')
                end
                
            end
            
            % CONTRIBUTION OF WEIGHTS
            %  USE_symmetric_form = 1;
            if USE_WEIGHT_TERM_k == 1
                %     if  USE_symmetric_form == 1
                %         K(DOFl, DOFl) =  K(DOFl, DOFl) + 0.5*(Kred_w_LL+Kred_w_LL') ;
                %   else
                
                K(DOFl, DOFl) =  K(DOFl, DOFl) + Kred_w_LL ;
                %  end
            end
            
            if CHECK_convergence == 1
                warning('Borrar esto ! ')
                eigAFTERw = eig(K(DOFl,DOFl))
                
                if any(eigAFTERw<0)
                    disp('Negative eigenvalue after WEIGHTS contribution')
                end
                
            end
            
            
            if DATA.ISDYNAMIC == 1
                % Dynamic part of the Jacobian
                K  =  KstifflDynamicPart(VAR,OPERFE,DATA.TIME_INTEGRATION_PARAMETERS,DATA.DeltaT,K) ;
            end
            
            if ~isempty(OPERFE.HYDRO)
                % Contribution of hydrostatic forces
                % \FextPR{}{}   = - \NbST{T}  \wSTft{}
                % diag(kPRESSst)*Lbool
                % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/DynamicFloatingBodies.tex
                % ------------------------
                K = K - OPERFE.HYDRO.NbST_w'*(ConvertBlockDiag_general(kPRESSst,DATA.MESH.ndim,OPERFE.HYDRO.irowsNDIM,OPERFE.HYDRO.icolsNDIM)*OPERFE.HYDRO.Lbool) ;
            end
            
            if isempty(OPERFE.DOFm)
                if DATA.SOLVER_IS_ITERATIVE == 1
                    tolCONJG =[] ;
                    niterCONJG =[] ;
                    delta_dL=    pcg(K(DOFl,DOFl),-VAR.RESID(DOFl),[],[]) ;
                else
                    delta_dL = -K(DOFl,DOFl)\VAR.RESID(DOFl);
                    
                    %                     disp('borrar esto')
                    %                     opts = struct('tol', 1e-10, 'maxit', 10000);
                    %                     ddd = eigs(K(DOFl,DOFl), 10, 'smallestabs', opts);
                    %                     disp(['Smallest EIG = ',num2str(ddd(1))])
                    %                     if ddd(1)<= 0
                    %                         disp('BUCKLING')
                    %                         ddd
                    %                     end
                    %
                    %
                    % %                     [EEE,VVV] =  eig(full(K(DOFl,DOFl)));
                    % %                     VVV = diag(VVV) ;
                    % %                     VVV_smallest =VVV(end);
                    % %                     disp('---------------------------')
                    %
                    
                end
                VAR.DISP(DOFl) =  VAR.DISP(DOFl) + delta_dL ;
                
                % CHECK-IN for the DAMAGE model 
                % -----------------------------
             %   VAR = CheckInSignDamageModel(VAR,DATA) ; 
                
                
                normd = norm(delta_dL) ;
                
            else
                % prescribed boundary conditions
                % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/10_PeriodicQUADL.mlx
                if DATA.SOLVER_IS_ITERATIVE == 1
                    tolCONJG =[] ;
                    niterCONJG =[] ;
                    delta_dL=    pcg(OPERFE.A'*(K*OPERFE.A),-VAR.RESID(DOFl),[],[]) ;
                else
                    delta_dL = -(OPERFE.A'*(K*OPERFE.A))\(OPERFE.A'*VAR.RESID);
                end
                VAR.DISP(DOFl) =  VAR.DISP(DOFl) + delta_dL ;
                VAR.DISP(OPERFE.DOFr) =  OPERFE.G*VAR.DISP(OPERFE.DOFm) + OPERFE.uBAR ;
                
                
                normd = norm(delta_dL) ;
                
                
                
            end
            
            
            
        end
        kiter = kiter + 1;
        
        %         if kiter == DATA.NEWTON_RAPHSON.NMAXiter && GEO_CONTRIBUTION == 1
        %             GEO_CONTRIBUTION = 0 ;
        %             DATA.NEWTON_RAPHSON.NMAXiter = 2*DATA.NEWTON_RAPHSON.NMAXiter;
        %             disp('Proceeding without geometric contribution')
        %         end
        
        %         if kiter > DATA.NEWTON_RAPHSON.NMAXiter && ~isempty(OPERFE.wSTs_cluster) &&  iterCHANGE_ECM_rule <=nmax_iterCHANGE_ECM_rule
        %             q = VAR.DISP(OPERFE.wSTs_cluster.IndexDOFl_q);
        %             [~, idx] = min(abs(OPERFE.wSTs_cluster.q - q));
        %             OPERFE.wSTs = OPERFE.wSTs_cluster.Values(:,idx) ;
        %             disp('---------------')
        %             disp(['ITER_IN_CHANGE_ECM = ',num2str(iterCHANGE_ECM_rule)])
        %             disp(['  Local indexes points nonzero weights =',num2str(find(OPERFE.wSTs~=0)')])
        %             disp('---------------')
        %             kiter = 1;
        %             iterCHANGE_ECM_rule= iterCHANGE_ECM_rule+1;
        %
        %         end
        
        
    end
end