function [VAR,CONVERGED,DATA,dQrot] = NewtonRapshonStaticLarge_COROT2(DATA,OPERFE,VAR,MATPRO,DOFl,KcLINq,BstQ,D_QrotALL,dQrot,D_QrotINIall)
% Adaptation of NewtonRapshonStaticLarge.m, corotational approach, v2
% See  /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTexAPPV.tex
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/04_UPDrotUNC.mlx
% JAHO, 02-DEC-2024, Honest Greens, Pedralbes, Barcelona.
% ----------------------------------------------------
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

% INITIALIZATION dQrot




while  kiter <=DATA.NEWTON_RAPHSON.NMAXiter
    
    % Update rotational displacements
    % ---------------------------------------------------------------------
    %---------------------------------------------------------------------------
    % Diagonalization increment of rotation
    D_dQrotALL = QrotINIall_GLOBAL(dQrot,OPERFE.INDEXsparseROTmat) ;    
    % Latex expression
    % \dCqLOC & =   \DiagC{\QrotALL^T}   \Par{\DiagC{ \dQrotALL  \QrotINTall} \XcALL - \XcALL}
    %    \\& = \DiagC{\QrotALL^T}   \Par{\DiagC{ \dQrotALL  \QrotALL \QrotINIall^T} \XcALL - \XcALL}
    dCq = D_dQrotALL * (D_QrotALL * (D_QrotINIall' * OPERFE.XcALL)) - OPERFE.XcALL ;  % Global reference system
    dCqLOC = D_QrotALL'*dCq ; % Local reference system
    
    
    % -------------------------------
    % Residual forces from displacements  (STATIC )
    % Before 9-Feb-2022
    %     [VAR.GLSTRAINS,VAR.PK2STRESS,celastST,VAR.RESID,Fint,FgradST,VAR.PK1STRESS,detFgrad] =...
    %         ResidualFromDisplacements(OPERFE,VAR.DISP,MATPRO,DATA,VAR.FEXT) ;%
    % After 9-Feb-2022
    DATA.kiter = kiter;
    [VAR,celastST,FgradST,detFgrad] =  ResidualFromDisplacementsVAR_COROT(OPERFE,VAR,MATPRO,DATA,VARint_n,...
        dCqLOC,KcLINq,BstQ,D_QrotALL) ;
    
    % RESIDUAL COMPATIBILITY CONDITION 
    % ---------------------------------
    % \ResCOMPall =    \DiagC{\YcmpD} \dClocINCRE
    % where \dClocINCRE \defeq  \DiagC{\QrotALL}^T  \LboolCall \dC -   \dCqLOC
    VAR.dClocINCRE = D_QrotALL'*(OPERFE.LboolCall*VAR.DISP)-dCqLOC ; 
    ResCOMPall = OPERFE.D_YcmpD*VAR.dClocINCRE ; 
     
    
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
        [~, ~,CONVERGED]=  CheckConvergenceLSTR_corot2(VAR_FINT_DOFL,VAR_FEXT_DOFl,VAR_RESID_DOFL,kiter,normd,DATA,ResCOMPall) ;
        
        
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
            
            K = GetStiffnessMatrix_COROT(OPERFE,DATA,VAR,FgradST,celastST,KcLINq,BstQ) ;
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
            Delta_d = zeros(size(VAR.DISP)) ; 
            Delta_d(DOFl) = delta_dL ; 
            dQrot =   RotUpdateCorot2(OPERFE,dQrot,DATA,ResCOMPall,Delta_d,D_QrotALL) ; 
            % ----------------------------------------------------------
            
            
            
        end
        kiter = kiter + 1;
    end
end



