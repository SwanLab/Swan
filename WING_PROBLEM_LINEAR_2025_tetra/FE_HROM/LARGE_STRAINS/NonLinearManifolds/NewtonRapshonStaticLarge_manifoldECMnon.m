function [VAR,CONVERGED,DATA] = NewtonRapshonStaticLarge_manifoldECMnon(DATA,OPERFE,VAR,MATPRO,DOFl)
%--------------------------------------------------------------------------
% function [VAR, CONVERGED, DATA] = ...
%     NewtonRapshonStaticLarge_manifoldECMnon(DATA, OPERFE, VAR, MATPRO, DOFl)
%
% PURPOSE:
%   Solves a static nonlinear equilibrium problem in reduced-order models 
%   where the solution lies on a nonlinear low-dimensional manifold τ(q).
%   The unknowns are reduced coordinates `q`, and Newton–Raphson iterations 
%   are performed in the tangent space of the manifold.
%
%   This function implements the tangent-manifold Newton solver described in 
%   Section 9.2 of *MLEARNstruct_1.pdf*, including internal force evaluation, 
%   residual projection, and reduced tangent assembly using ECM weights and 
%   extrapolated slave contributions.
%
% KEY FEATURES:
%   - Newton–Raphson iterations in reduced coordinates (q-space)
%   - Projection of residual and tangent matrix onto manifold tangent space
%   - Nonlinear extrapolation of internal forces via η_NON
%   - Optional use of geometric stiffness (via τ″(q)) for improved accuracy
%   - Automatic treatment of Dirichlet boundary conditions (prescribed DOFs)
%
% INPUTS:
% --------
%   DATA     : struct with model settings, tolerances, solver flags,
%              and internal variable names for history-dependent materials.
%
%   OPERFE   : struct with reduced-order operators and mappings:
%              - tauNON(q): nonlinear mapping τ(q)
%              - tauNONder(q): Jacobian ∂τ/∂q
%              - tauNONder2(q): second derivatives for geometric stiffness
%              - etaNON, etaNONder: nonlinear extrapolation of internal force
%              - DOFl, DOFr, DOFm, A, G, uBAR: DOF maps for BC handling
%
%   VAR      : struct with current solution:
%              - DISP: current displacements
%              - FEXT_extended: full-space external force
%              - other fields updated during iterations (stresses, residuals, etc.)
%
%   MATPRO   : struct with constitutive law and material properties
%
%   DOFl     : indices of free (independent) degrees of freedom
%
% OUTPUTS:
% ---------
%   VAR        : updated structure containing converged displacement and stresses
%   CONVERGED  : logical flag indicating if convergence was achieved
%   DATA       : updated data structure with iteration history
%
% ALGORITHM OVERVIEW:
% -------------------
%   1. Extract reduced coordinate q from VAR.DISP(DOFl)
%   2. Evaluate τ(q), τ′(q), and τ″(q) (if available)
%   3. Build full displacement field: u = [τ(q); qR]
%   4. Compute internal forces and residual via `ResidualFromDisplacementsVARnonECM`
%   5. Assemble projected tangent matrix using `GetStiffnessMatrix_nonECM`
%   6. Solve reduced linear system: Δq ← -K⁻¹·r
%   7. Update q and check convergence
%   8. On convergence: compute Von Mises stresses, strain energy, optional dynamics
%
% FORMULATION REFERENCES (MLEARNstruct_1.pdf):
%   - τ(q) mapping:               Eq. (176)
%   - Projected residual:        Eq. (182)
%   - Reduced tangent matrix:    Eq. (184)
%   - Equivalent master weights: Eq. (188)
%   - Geometric correction:      Eq. (191–193)
%
% REMARKS:
% ---------
%   - Compatible with linear or nonlinear materials, with optional dynamics
%   - Dirichlet BCs handled using constraint projection matrices (A, G)
%   - Uses Empirical Cubature Method (ECM) with master/slave integration
%   - Designed for static problems; dynamic extension flagged as TODO
%
% DEPENDENCIES:
%   - ResidualFromDisplacementsVARnonECM
%   - InternalForcesECMnon2
%   - GetStiffnessMatrix_nonECM
%   - CheckConvergenceLSTR
%   - VonMisesCauchyStresses
%   - KineticAndStrainEnergyGet
%
% AUTHOR:
%   Joaquín A. Hernández Ortega, Universitat Politècnica de Catalunya (UPC) – CIMNE
%   Last revised: 22 July 2025, Pedralbes (Barcelona)
%--------------------------------------------------------------------------

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
    tauNON_q= OPERFE.tauNON(qL);
    tauNONder_q = OPERFE.tauNONder(qL);
    if ~isempty( OPERFE.tauNONder2)
        tauNONder2_q = OPERFE.tauNONder2(qL);
        % this is only for the unconstrained DOFs
    else
        tauNONder2_q = [] ;
    end
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