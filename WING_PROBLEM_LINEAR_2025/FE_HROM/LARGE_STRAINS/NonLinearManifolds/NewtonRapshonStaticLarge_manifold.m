function [VAR,CONVERGED,DATA] = NewtonRapshonStaticLarge_manifold(DATA,OPERFE,VAR,MATPRO,DOFl)
%--------------------------------------------------------------------------
% function [VAR, CONVERGED, DATA] = NewtonRapshonStaticLarge_manifold(DATA, OPERFE, VAR, MATPRO, DOFl)
%
% PURPOSE:
%   Solves a static equilibrium problem on a nonlinear solution manifold
%   τ(q) using Newton–Raphson iterations projected onto the tangent space
%   of the manifold. This implements the tangent-manifold formulation
%   described in Section 9.2 of *MLEARNstruct_1.pdf*.
%
%   Instead of solving the full system in nodal coordinates, the unknowns
%   are generalized coordinates q that parametrize the solution manifold.
%   At each iteration, the residual is projected onto the tangent space
%   of the manifold using the Jacobian τ′(q), and the update Δq is computed
%   from the reduced system. Convergence is checked in the reduced space.
%
% INPUTS:
%   - DATA     : structure with simulation controls, tolerances, time settings
%   - OPERFE   : operators and mappings, including:
%                   * tauNON(q):      nonlinear map τ(q)
%                   * tauNONder(q):   tangent Jacobian τ′(q)
%                   * tauNONder2(q):  second derivative for geometric stiffness
%                   * DOFr, DOFm:     index sets for prescribed Dirichlet constraints
%   - VAR      : structure with displacement, forces, stresses, internal vars
%   - MATPRO   : constitutive model parameters
%   - DOFl     : logical or index array of free degrees of freedom
%
% OUTPUTS:
%   - VAR        : updated solution structure after convergence
%   - CONVERGED  : flag indicating whether Newton iteration converged
%   - DATA       : updated DATA (e.g., iteration history, residual norms)
%
% ALGORITHM OVERVIEW:
%   1. Extract reduced coordinates q from VAR.DISP(DOFl)
%   2. Evaluate τ(q), τ′(q), and τ″(q) if needed
%   3. Construct full displacement field u = τ(q) + prescribed terms
%   4. Compute residual R(u) and internal forces using constitutive model
%   5. Project residual and tangent matrix to reduced space via τ′(q)ᵗ
%   6. Solve reduced Newton system and update q ← q + Δq
%   7. Check convergence (in reduced residual or update norm)
%   8. If converged, post-process Von Mises stress and energies
%
% FORMULATION:
%   - Full-space residual:             r = f_ext - f_int(u)
%   - Projected residual (Eq. 182):    r_ROM = τ′(q)ᵗ · r
%   - Projected tangent (Eq. 184):     K_ROM = τ′(q)ᵗ · K(u) · τ′(q)
%   - Geometric correction:            ∑_{i} τ″_i(q) · r_i
%
% REFERENCES:
%   *MLEARNstruct_1.pdf*, Section 9.2, pp. 71–77
%     - Eq. (176): definition of τ(q)
%     - Eq. (182): residual projection
%     - Eq. (184): tangent matrix projection
%     - Eq. (185)–(193): generalization to include hydrostatics, dynamics, etc.
%
% ADDITIONAL FEATURES:
%   - prescribed BCs handled via partitioned DOF sets and matrices A, G, uBAR
%   - Optional geometric stiffness via τ″(q) when second derivatives are available
%   - Can be extended to dynamic or pressure-driven systems
%
% DEPENDENCIES:
%   - ResidualFromDisplacementsVAR
%   - GetStiffnessMatrix
%   - CheckConvergenceLSTR
%   - VonMisesCauchyStresses
%   - UpdateVelocityAcceleration
%   - KineticAndStrainEnergyGet
%   - PressureForcesFromDisplacements (optional)
%   - ResidualDynamicPart (optional)
%
% AUTHOR:
%   Joaquín A. Hernández Ortega, UPC – CIMNE
%   Last updated: 2 July 2025, Honest Greens – Pedralbes, Barcelona
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