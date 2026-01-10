function [VAR,CONVERGED,DATA,Qrot] = NewtonRapshonStaticLarge_COROT_LRss(DATA,OPERFE,VAR,MATPRO,DOFl,Qrot,Delta_dR,DOFr)
% ---------------------------------------------------------------------------------------------------
% FUNCTION: NewtonRapshonStaticLarge_COROT_LRss
% ---------------------------------------------------------------------------------------------------
% PURPOSE:
%   Newton–Raphson solver tailored for geometrically nonlinear problems in the **small strain, large rotation**
%   (SSLR) or **large strain, small rotation** (LSSR) regime using a **corotational framework**.
%   This function updates both displacements and rotation matrices to enforce equilibrium in the local frame.
%
% USAGE:
%   [VAR, CONVERGED, DATA, Qrot] = NewtonRapshonStaticLarge_COROT_LRss(DATA, OPERFE, VAR, MATPRO, DOFl, Qrot, Delta_dR, DOFr)
%
% INPUTS:
%   - DATA       : Simulation control structure with settings for Newton solver, convergence, dynamic options, etc.
%   - OPERFE     : Structure containing finite element operators, geometric mappings, rotation mappings, etc.
%   - VAR        : Structure with evolving field variables: displacements, internal variables, residuals, stresses.
%   - MATPRO     : Material property structure.
%   - DOFl       : Indices of free (active) degrees of freedom for displacements.
%   - Qrot       : Structure (or matrix) holding the rotation matrices per element or node (diagonal blocks).
%   - Delta_dR   : Prescribed macroscopic displacement increment (used for rotational update).
%   - DOFr       : Indices of constrained DOFs (used in affine BC application).
%
% OUTPUTS:
%   - VAR        : Updated field variable structure (displacements, stresses, forces, energies).
%   - CONVERGED  : Flag (1 if Newton converged, 0 otherwise).
%   - DATA       : Updated data structure (iteration tracking, snapshot storage, etc.).
%   - Qrot       : Updated rotation matrices for all nodes/elements after convergence.
%
% FUNCTIONALITY:
%   - Builds the rotated residual using a corotated frame (local/global projection via Qrot).
%   - Assembles tangent stiffness matrix in local rotated coordinates, including geometric and coupling terms.
%   - Updates nodal rotations using rotation update law (`RotUpdateCorot_LR`).
%   - Evaluates convergence using residual and force norm criteria.
%   - Supports:
%       * Hydrostatic forces (`OPERFE.HYDRO`)
%       * Dynamic residual and tangent contributions (`ISDYNAMIC == 1`)
%       * Affine or periodic boundary conditions (via mapping operators A and G)
%
% SPECIAL FEATURES:
%   - Designed for mixed kinematics: small strains with large rigid-body-like rotations (e.g., beams, shells, micropolar continua).
%   - Integrates with EIFEM framework and supports EIFEM-compatible rotation handling.
%   - Modular handling of rotation–displacement coupling via matrix operators (`LboolCallQ`, `PdownsRBcoupROTc_i`).
%
% REFERENCES:
%   - See derivations and applications in:
%     * /TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
%     * /PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf
%
% AUTHOR:
%   Joaquín A. Hernández, UPC – CIMNE
%   Date: 7-Feb-2025, Friday
%   Comments by ChatGPT-4, 12-May-2025
%
% DEPENDENCIES:
%   - ResidualStaticCOROTfinal
%   - GetStiffnessMatrix_COROT_LRss
%   - VonMisesCauchyStresses_COROT
%   - KineticAndStrainEnergyGet_COROT2
%   - RotUpdateCorot_LR
%   - UpdateVelocityAcceleration
%   - CheckConvergenceLSTR
%   - DefaultField (utility)
%
% ---------------------------------------------------------------------------------------------------

% Adaptation of NewtonRapshonStaticLarge.m,  
% See  Small strains/Large rotations (or Small rotations/Large strains)
% JAHO, 7-feb-2025, FRIDAY, 14:31,   Barcelona.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
% and
%/home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf
% ----------------------------------------------------
if nargin == 0
    load('tmp3.mat')
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
    
    
    % Diagonal global rotation matrix:  \DiagC{\QrotALL}:  
    D_QrotALL = QrotINIall_GLOBAL(Qrot,OPERFE.INDEXsparseROTmat) ;    
    %  Rotated Boolean operator:  \LboolCallQ =   \DiagC{\QrotALL}^T \LboolCall 
    LboolCallQ = D_QrotALL'*OPERFE.LboolCall; 
    % Rotational displacements in local coordinates 
    % \dCqLOC  =  \XcALLloc - \DiagC{\QrotALL}^T \XcALL
    dCqLOC = OPERFE.XcALLloc - D_QrotALL'*OPERFE.XcALL ;     
    % Static residual 
     [VAR,celastST,FgradST,detFgrad,KcGEOunassemGLOloc,PdownsRBcoupROTc_i] = ...
         ResidualStaticCOROTfinal(OPERFE,D_QrotALL,VAR,LboolCallQ,MATPRO,VARint_n,DATA,dCqLOC,Qrot) ; 
    
    if DATA.ISDYNAMIC == 1 && ~isempty(VAR.RESID)
        % Dynamic part of the residual
        VAR =  ResidualDynamicPart(VAR,OPERFE,DATA.TIME_INTEGRATION_PARAMETERS,DATA.DeltaT) ;
    end
    if ~isempty(OPERFE.HYDRO)
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
        disp('Not converged ..., because of Negative JAcobians')
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
                VonMisesCauchyStresses_COROT(PK1STRESS_forCauchy,FgradST,DATA.MESH.ndim,detFgrad,DATA) ;
            if DATA.ISDYNAMIC == 1
                % Update velocities and accelerations
                VAR =  UpdateVelocityAcceleration(VAR,DATA.TIME_INTEGRATION_PARAMETERS,DATA.DeltaT) ;
                %  VAR = KineticAndStrainEnergyGet(VAR,OPERFE,DATA,MATPRO) ;
                %   6-May-2023
                %
            end
            VAR = KineticAndStrainEnergyGet_COROT2(VAR,OPERFE,DATA,MATPRO,FgradST,dCqLOC,D_QrotALL) ;
            break
        else
            %  Updating d_l by solving  Delta_d_L =  -K(DOFl,DOFl)\Res(DOFl)
            % 8. Assembly Stiffness Matrix (static part)
            
            K = GetStiffnessMatrix_COROT_LRss(OPERFE,DATA,VAR,FgradST,celastST,LboolCallQ,D_QrotALL,KcGEOunassemGLOloc) ;
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
                    
%                                         disp('borrar esto')
%                     opts = struct('tol', 1e-10, 'maxit', 10000);
%                     ddd = eigs(K(DOFl,DOFl), 10, 'smallestabs', opts);
%                     disp(['Smallest EIG = ',num2str(ddd(1))])
%                     if ddd(1)<= 0
%                         disp('BUCKLING')
%                         ddd
%                     end
                    
                    
                end
                VAR.DISP(DOFl) =  VAR.DISP(DOFl) + delta_dL ;
                normd = norm(delta_dL) ;
                
            else
                % Affine boundary conditions
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
            
            % ROTATIONS UPDATE ------------------------------------------
            
            
             Qrot =   RotUpdateCorot_LR(OPERFE,Qrot,DATA,delta_dL,LboolCallQ,Delta_dR,kiter,VAR,DOFl,DOFr,PdownsRBcoupROTc_i) ; 




          %  [Qrot,norm_angle_rads] = UpdateRotationsEIFEM2(OPERFE,Qrot,VAR,DATA,iterROTATIONS,dQrot) ;
            % ----------------------------------------------------------
            
            
            
        end
        kiter = kiter + 1;
    end
end



