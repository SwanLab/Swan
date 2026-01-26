function [VAR,CONVERGED,DATA] = NewtonRapshonStaticLarge_manifoldFST(DATA,OPERFE,VAR,MATPRO,DOFl)
%--------------------------------------------------------------------------
% NewtonRapshonStaticLarge_manifoldFST
%
% PURPOSE:
%   Faster, updated variant of `NewtonRapshonStaticLarge_manifold` for
%   solving static equilibrium on a nonlinear manifold τ(q) via
%   Newton–Raphson in the manifold’s tangent space. This FST version
%   changes HOW modal amplitudes are obtained from the generalized
%   coordinates q: it uses a precompiled evaluator (no anonymous functions),
%   enabling vectorized, low-overhead evaluations of τ(q), τ′(q), τ″(q).
%
% WHAT CHANGED VS ORIGINAL:
%   • Replaces anonymous-function evaluators for τ, τ′, τ″ with a
%     precompiled structure (e.g., `DATA.DATA_evaluateTAU_and_DER`) created
%     during the FAST offline stage. This preserves theory but accelerates
%     per-iteration evaluations.
%
% METHOD OVERVIEW (unchanged theory):
%   Given the ROM parametrization u = τ(q) with tangent τ′(q),
%   Newton iterations are carried out in reduced space:
%     1) Form full residual r = f_ext − f_int(u) at current u = τ(q).
%     2) Project to reduced residual  r_ROM = τ′(q)ᵀ r.
%     3) Project tangent            K_ROM = τ′(q)ᵀ K(u) τ′(q).
%     4) (Optional) add geometric term with τ″(q)·r (curvature contribution).
%     5) Solve K_ROM Δq = −r_ROM and update q ← q + Δq.
%     6) Check convergence in reduced space; postprocess stresses/energies.
%
% INPUTS:
%   - DATA     : Simulation controls, tolerances, flags (includes pointer to
%                the τ/derivatives evaluator struct).
%   - OPERFE   : Operators & DOF sets (DOFl/DOFr/DOFm), mass/stress ops,
%                optional hydro, etc.
%   - VAR      : Current state (DISP, RESID, FINT/FEXT, stresses, etc.).
%   - MATPRO   : Constitutive model parameters.
%   - DOFl     : Index/logic of free reduced DOFs.
%
% OUTPUTS:
%   - VAR       : Updated state after convergence (stresses, energies, etc.).
%   - CONVERGED : Logical flag (1 if converged, 0 otherwise).
%   - DATA      : Updated with iteration history / norms / snapshots.
%
% FEATURES:
%   - Enforces affine/periodic BCs via partition (DOFl/DOFr/DOFm) and maps A,G.
%   - Optional geometric stiffness contribution via τ″(q).
%   - Compatible with small- or large-strain formulations (as configured).
%   - Snapshot collection each step for later reconstruction/analysis.
%
% REFERENCES (theory unchanged from original):
%   * MLEARNstruct_1.pdf, §9.2:
%       - Eq. (176) τ(q) definition
%       - Eq. (182) residual projection τ′(q)ᵀ r
%       - Eq. (184) tangent projection τ′(q)ᵀ K τ′(q)
%       - Eq. (185)–(193) geometric/dynamic extensions
%
% DEPENDENCIES:
%   - ResidualFromDisplacementsVAR
%   - GetStiffnessMatrix
%   - CheckConvergenceLSTR
%   - VonMisesCauchyStresses
%   - UpdateVelocityAcceleration (if dynamic extension used)
%   - KineticAndStrainEnergyGet
%
% HISTORY (CREATION / UPDATES):
%   • Created  : 02-Jul-2025, Honest Greens – Pedralbes (Barcelona)
%                (as `NewtonRapshonStaticLarge_manifold`)
%   • Updated  : 11-Aug-2025, Molinos Marfagones (Cartagena)
%                - FST variant: uses precompiled τ/derivatives evaluator
%                  (no anonymous functions); same formulation, faster evals.
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


while  kiter <=DATA.NEWTON_RAPHSON.NMAXiter
    
    if  ~isempty(OPERFE.wSTs_cluster)
    q = VAR.DISP(OPERFE.wSTs_cluster.IndexDOFl_q); 
w  = feval(OPERFE.wSTs_cluster.DATA_regress.nameFunctionEvaluate,q,OPERFE.wSTs_cluster.DATA_regress) ;  
OPERFE.wSTs = w; 
 
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
    [VAR,celastST,FgradST,detFgrad] =  ResidualFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
    % Now we have to project some of these variables on the tangent manifold
    % To this end, we need tauNONallDER_q
    % This is a matrix of size nmodes_extended x ndof
    tauNONallDER_q_RR = speye(length(OPERFE.DOFr),length(OPERFE.DOFr)) ;
    tauNONallDER_q_LR = sparse(size(tauNONder_q,1),length(OPERFE.DOFr));
    tauNONallDER_q_RL =  sparse(length(OPERFE.DOFr),size(tauNONder_q,2));
    
    tauNONallDER_q = [tauNONder_q,tauNONallDER_q_LR
        tauNONallDER_q_RL, tauNONallDER_q_RR ] ;
    
    % Projection onto the reduced coefficients
    VAR.RESID_extend = VAR.RESID ;
    VAR.FINT = tauNONallDER_q'*VAR.FINT ;
    VAR.FEXT = tauNONallDER_q'*VAR.FEXT ;
    VAR.RESID = tauNONallDER_q'*VAR.RESID ;
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
            K = GetStiffnessMatrix(OPERFE, DATA, VAR, FgradST, celastST) ;
            
            %--------------------------------------------------------------------------
            % 2. Project the full stiffness matrix onto the tangent space of the nonlinear manifold
            %    Following Eq. (184):   K_ROM = τ′(q)^T · K · τ′(q)
            %    This yields the reduced tangent matrix governing the update of q in the ROM.
            K = tauNONallDER_q' * K * tauNONallDER_q ;
            
            %--------------------------------------------------------------------------
            % 3. Optional: Add second-order geometric correction to the stiffness matrix
            %    If τ″(q) is available, use Eq. (185) from the reference to include
            %    geometric stiffness arising from the curvature of the manifold.
            %    This term is: ∑_i [ τ″_i(q) * r_i(u) ]  (scalar-weighted Hessian contributions)
          %  tauNONder2_q = [] ; 
           % disp('borrar')
           if ~isempty(tauNONder2_q)
                K_geo = 0 ;
                
                % Loop over generalized DOFs (excluding prescribed DOFs)
                for imodesEXT = 1:(length(VAR.RESID_extend) - length(OPERFE.DOFr))
                    % Geometric stiffness contribution: second derivative * residual component
                    K_geo = K_geo + tauNONder2_q(imodesEXT) * VAR.RESID_extend(imodesEXT);
                end
                
                % Add geometric stiffness contribution to the reduced tangent matrix
                % Only affecting the reduced DOFs (DOFl)
                K(DOFl, DOFl) = K(DOFl, DOFl) + K_geo ;
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