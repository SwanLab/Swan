function [VAR,CONVERGED,DATA] = NewtonRapshonStaticLarge_manifoldECMnonFST(DATA,OPERFE,VAR,MATPRO,DOFl)
%--------------------------------------------------------------------------
% NewtonRapshonStaticLarge_manifoldECMnonFST
%
% PURPOSE:
%   Faster, updated variant of `NewtonRapshonStaticLarge_manifoldECMnon`
%   for static manifold-ROM solves with ECM hyperreduction (master/slave).
%   This FST version changes HOW η (master→slave reconstruction of internal
%   work densities) and its derivatives are evaluated:
%     • Replaces anonymous-function evaluators with a precompiled,
%       vectorized evaluator structure saved offline
%       (e.g., `OPERFE.DATA_regress_eta_der` produced by
%        `DiscreteECM_adaptWEIGHTSfst` → `BsplinesLeastSquares_fastECM`).
%     • Supports equivalent master weights `wMSTeq` and reduced residuals
%       `ResREDeqMST` returned by the residual routine for a leaner Newton step.
%
% METHOD OVERVIEW (theory unchanged from original):
%   1) Use reduced coordinates q (free DOFs) and construct u = [τ(q); qR].
%   2) Compute residual and internal forces on the ECM points, with slave
%      contributions reconstructed via η(q_mst).
%   3) Project residual/tangent to reduced space; assemble K using the ECM,
%      equivalent weights, and (optionally) geometric terms from τ″(q).
%   4) Solve K Δq = −r, update q, check convergence.
%   5) On convergence, post-process stresses/energies.
%
% WHAT CHANGED VS ORIGINAL (`…ECMnon`):
%   • η and its derivatives are evaluated through a precompiled regression
%     object (`DATA_regress_eta_der`) → eliminates per-iteration anonymous
%     functions and enables vectorized evaluation.
%   • Integrates equivalent master-point weights `wMSTeq` and reduced
%     equations `ResREDeqMST` in stiffness/residual assembly.
%
% INPUTS:
%   - DATA     : Simulation controls and options (tolerances, flags, step data).
%   - OPERFE   : Reduced operators and mappings:
%                   • DOFl/DOFr(/DOFm), A, G, uBAR
%                   • Bst, wSTs, (optional) HYDRO
%                   • DATA_regress_eta_der  ← precompiled η/η′/η″ evaluator
%   - VAR      : Current state (DISP, RESID, FINT/FEXT, stresses, etc.).
%   - MATPRO   : Material/constitutive parameters.
%   - DOFl     : Indices (or logical mask) of free reduced DOFs.
%
% OUTPUTS:
%   - VAR        : Updated state after convergence (stresses, energies, etc.).
%   - CONVERGED  : Logical flag (1 if converged, 0 otherwise).
%   - DATA       : Updated with iteration history / snapshots.
%
% FORMULATION REFERENCES (MLEARNstruct_1.pdf):
%   • §9.2  – Manifold ROM: τ(q), residual/tangent projections.
%   • §9.3  – ECM hyperreduction & master/slave reconstruction (η).
%   • Eqs. (176), (182), (184), (188), (191–193).
%
% DEPENDENCIES:
%   - ResidualFromDisplacementsVARnonECM
%   - GetStiffnessMatrix_nonECM
%   - CheckConvergenceLSTR
%   - VonMisesCauchyStresses
%   - KineticAndStrainEnergyGet
%   - (optional) UpdateVelocityAcceleration, PressureForcesFromDisplacements
%
% NOTES:
%   - Current implementation assumes 1 master integration point in stiffness
%     assembly (extend as needed).
%   - Optional numerical derivative path available via DATA flags for
%     debugging/verification.
%   - Compatible with affine/periodic BCs via A, G, DOFm.
%
% HISTORY (CREATION / UPDATES):
%   • Created  : 22-Jul-2025 (as `NewtonRapshonStaticLarge_manifoldECMnon`)
%   • Updated  : 12-Aug-2025, Molinos Marfagones (Cartagena)
%                - FST variant: precompiled η/η′/η″ evaluator, equivalent
%                  master weights, streamlined Newton assembly.
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC/CIMNE)
%   Comments updated by ChatGPT, 12-Aug-2025
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

%

if nargin == 0
    load('tmp1.mat')
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


DATA.EvaluateDerivativeNumerically = 0;
DATA.IntervalEvaluationDerivative = 1e-6;


while  kiter <=DATA.NEWTON_RAPHSON.NMAXiter
    
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
    % Now we create the displacement using "extended" coordinate
    VAR.DISP = [tauNON_q;qR] ;
    % Now we can invoke the function that computes the residual
    DATA.kiter = kiter;
    VAR.FEXT = VAR.FEXT_extended;
    % Now we have to project some of these variables on the tangent manifold
    % To this end, we need tauNONallDER_q
    % This is a matrix of size nmodes_extended x ndof
    tauNONallDER_q_RR = speye(length(OPERFE.DOFr),length(OPERFE.DOFr)) ;
    tauNONallDER_q_LR = sparse(size(tauNONder_q,1),length(OPERFE.DOFr));
    tauNONallDER_q_RL =  sparse(length(OPERFE.DOFr),size(tauNONder_q,2));
    tauNONallDER_q = [tauNONder_q,tauNONallDER_q_LR
        tauNONallDER_q_RL, tauNONallDER_q_RR ] ;    % The function below is for evaluating internal forces and residual
    
    [VAR,celastST,FgradST,detFgrad,fINTredNON_mst_DOFl,wMSTeq,ResREDeqMST] =...
        ResidualFromDisplacementsVARnonECM(OPERFE,VAR,MATPRO,DATA,VARint_n,tauNONallDER_q) ;
    
    VAR.EquivalentWeightsECM =  wMSTeq;
    
    
    if  DATA.EvaluateDerivativeNumerically >= 1;
        der_Fint_q = EvaluateResidualDerivative(DATA,OPERFE,qR,VAR,qL,MATPRO,VARint_n) ;
    else
        der_Fint_q = [] ;
    end
    
    % Projection onto the reduced coefficients
    % VAR.RESID_extend = VAR.RESID ;  Not used in this version
    %  VAR.FINT = tauNONallDER_q'*VAR.FINT ;
    
    % Go back to the reduced coordinates
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
        [~, ~,CONVERGED]=  CheckConvergenceLSTR(VAR_FINT_DOFL,VAR_FEXT_DOFl,VAR_RESID_DOFL,kiter,normd,DATA) ;
        
        
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
            VAR = KineticAndStrainEnergyGet(VAR,OPERFE,DATA,FgradST) ;  % 6-May-2023
            break
        else
            
            %--------------------------------------------------------------------------
            % 1. Compute full-space tangent stiffness matrix based on current deformation state
            %    This includes material and geometric contributions in the full displacement space.
            
            % On July 21st 2025, the only option available is when the
            % number of master integration points is equal to 1
            %
            mMST = length(OPERFE.wSTs) ;
            
            if  mMST == 1
                %  Bst_non = OPERFE.Bst*tauNONallDER_q ;
                
                K = GetStiffnessMatrix_nonECM(OPERFE, DATA, VAR, FgradST, celastST,wMSTeq,ResREDeqMST,tauNONallDER_q,tauNONder2_q) ;
                
                if DATA.EvaluateDerivativeNumerically == 1
                    K(1,1) = der_Fint_q ;
                end
                
            else
                error('Option not implemented yet (22-July-2025)')
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
    end
end