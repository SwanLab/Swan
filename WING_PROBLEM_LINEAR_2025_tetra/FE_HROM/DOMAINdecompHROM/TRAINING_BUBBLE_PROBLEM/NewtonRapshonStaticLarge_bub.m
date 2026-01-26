function [VAR,CONVERGED,DATA] = NewtonRapshonStaticLarge_bub(DATA,OPERFE,VAR,MATPRO,DISP_CONDITIONS)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/08_IndividualTRAINbub.mlx
if nargin == 0
    load('tmp3.mat')
end
%

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
    DATA.kiter = kiter;
    [VAR,celastST,FgradST,detFgrad] =  ResidualFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
    
%     if DATA.ISDYNAMIC == 1 && ~isempty(VAR.RESID)
%         % Dynamic part of the residual
%         VAR =  ResidualDynamicPart(VAR,OPERFE,DATA.TIME_INTEGRATION_PARAMETERS,DATA.DeltaT) ;
%     end
%     if ~isempty(OPERFE.HYDRO)
%         % Hydrostatic forces
%         % \FextPR{}{}   = - \NbST{T}  \wSTft{}  \tPRESSst}
%         % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/DynamicFloatingBodies.tex
%         % ------------------------
%         [VAR.FEXTpress,VAR.tPRESS,kPRESSst ]= PressureForcesFromDisplacements(DATA,OPERFE,VAR) ;
%         VAR.RESID = VAR.RESID - VAR.FEXTpress ;
%         VAR_FEXT_DOFl = VAR.FEXT(DOFl) + VAR.FEXTpress(DOFl);
%     else
%         VAR_FEXT_DOFl = VAR.FEXT(DOFl) ;
%     end
    
    % Update snapshots of iterative results
    DATA = UpdateIterativeSnapshot(DATA,VAR,1+kiter) ;
    
    if isempty(celastST) && DATA.INTERNAL_FORCES_USING_precomputed_Kstiff == 0
        disp('Not converged ..., because of Negative JAcobians')
        DATA = ReshapeIterativeSnapshot(DATA,VAR,kiter) ;
        break
    else
        % 7. Check convergence
        VAR_FINT_DOFl = DISP_CONDITIONS.A'*VAR.FINT; 
        VAR_FEXT_DOFl = zeros(size(VAR_FINT_DOFl)) ; 
        VAR_RESID_DOFl = VAR_FINT_DOFl ; 
        
        [criter_f, norm_res,CONVERGED]=  CheckConvergenceLSTR(VAR_FINT_DOFl,VAR_FEXT_DOFl,VAR_RESID_DOFl,kiter,normd,DATA) ;
     %   [criter_f, norm_res,CONVERGED]=  CheckConvergenceLSTR(VAR.FINT(DOFl),VAR_FEXT_DOFl,VAR.RESID(DOFl),kiter,normd,DATA) ;
        if CONVERGED == 1
            % Compute Von Mises Stress
            if DATA.SMALL_STRAIN_KINEMATICS == 0
                PK1STRESS_forCauchy  = VAR.PK1STRESS ;
            else
                PK1STRESS_forCauchy  = VAR.PK2STRESS ;
            end
            [ VAR.VONMISES_CAUCHY_STRESS,VAR.CAUCHY_STRESS ]= ...
                VonMisesCauchyStresses(PK1STRESS_forCauchy,FgradST,DATA.MESH.ndim,detFgrad,DATA) ;
%             if DATA.ISDYNAMIC == 1
%                 % Update velocities and accelerations
%                 VAR =  UpdateVelocityAcceleration(VAR,DATA.TIME_INTEGRATION_PARAMETERS,DATA.DeltaT) ;
%                %  VAR = KineticAndStrainEnergyGet(VAR,OPERFE,DATA,MATPRO) ;
%                %   6-May-2023
%                %  
%             end
            VAR = KineticAndStrainEnergyGet(VAR,OPERFE,DATA,MATPRO) ;  % 6-May-2023
            break
        else
            %  Updating d_l by solving  Delta_d_L =  -K(DOFl,DOFl)\Res(DOFl)
            % 8. Assembly Stiffness Matrix (static part)
            
            K = GetStiffnessMatrix_bub(OPERFE,DATA,VAR,FgradST,celastST) ;

            if DATA.SOLVER_IS_ITERATIVE == 1
%                 tolCONJG =[] ;
%                 niterCONJG =[] ;
                delta_BUB_dL=    pcg(K,-DISP_CONDITIONS.A'*VAR.RESID,[],[]) ;
            else
                delta_BUB_dL = -K\(DISP_CONDITIONS.A'*VAR.RESID);
            end
            VAR.BUB_DISP(DISP_CONDITIONS.DOFl) =  VAR.BUB_DISP(DISP_CONDITIONS.DOFl) + delta_BUB_dL ; 
            VAR.BUB_DISP(DISP_CONDITIONS.DOFr) = DISP_CONDITIONS.J*VAR.BUB_DISP(DISP_CONDITIONS.DOFm)  ; 
            
            VAR.DISP = VAR.ELAST_DISP    + VAR.BUB_DISP ;
            normd = norm(delta_BUB_dL) ;
            
        end
        kiter = kiter + 1;
    end
end