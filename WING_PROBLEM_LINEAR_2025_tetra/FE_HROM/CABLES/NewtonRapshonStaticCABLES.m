function [VAR,CONVERGED,DATA] = NewtonRapshonStaticCABLES(DATA,OPERFE,VAR,MATPRO,DOFl)
if nargin == 0
    load('tmp.mat')
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
    DATA.kiter = kiter;
    [VAR,celastST,FgradST,detFgrad] =  ResidualFromDisplacementsCABLES(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
    
   
    
    if DATA.ISDYNAMIC == 1 && ~isempty(VAR.RESID)
        % Dynamic part of the residual
        VAR =  ResidualDynamicPart(VAR,OPERFE,DATA.TIME_INTEGRATION_PARAMETERS,DATA.DeltaT) ;
        
       
        
    end
%     if ~isempty(OPERFE.HYDRO)
%         % Hydrostatic forces
%         % \FextPR{}{}   = - \NbST{T}  \wSTft{}  \tPRESSst}
%         % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/DynamicFloatingBodies.tex
%         % ------------------------
%         [VAR.FEXTpress,VAR.tPRESS,kPRESSst ]= PressureForcesFromDisplacements(DATA,OPERFE,VAR) ;
%         VAR.RESID = VAR.RESID - VAR.FEXTpress ;
%         VAR_FEXT_DOFl = VAR.FEXT(DOFl) + VAR.FEXTpress(DOFl);
%     else
        VAR_FEXT_DOFl = VAR.FEXT(DOFl) ;
 %   end
    
    % Update snapshots of iterative results
    DATA = UpdateIterativeSnapshot(DATA,VAR,1+kiter) ;
    
    if (isempty(celastST) ) 
       
        disp('Not converged ..., because of Negative JAcobians')
        DATA = ReshapeIterativeSnapshot(DATA,VAR,kiter) ;
        break
    else
        % 7. Check convergence
        [criter_f, norm_res,CONVERGED]=  CheckConvergenceLSTR(VAR.FINT(DOFl),VAR_FEXT_DOFl,VAR.RESID(DOFl),kiter,normd,DATA) ;
        if CONVERGED == 1
            % Compute Von Mises Stress
%             if DATA.SMALL_STRAIN_KINEMATICS == 0
%                 PK1STRESS_forCauchy  = VAR.PK1STRESS ;
%             else
%                 PK1STRESS_forCauchy  = VAR.PK2STRESS ;
%             end
%             [ VAR.VONMISES_CAUCHY_STRESS,VAR.CAUCHY_STRESS ]= ...
%                 VonMisesCauchyStresses(PK1STRESS_forCauchy,FgradST,DATA.MESH.ndim,detFgrad,DATA) ;
           ElementsNegative = find(VAR.TENSION <0) ; 
           if ~isempty(ElementsNegative)
               disp(['Points in compression = ',num2str(ElementsNegative')])
             %  warning('Converged solution with elements in compression')
           end

 
            if DATA.ISDYNAMIC == 1
                % Update velocities and accelerations
                VAR =  UpdateVelocityAcceleration(VAR,DATA.TIME_INTEGRATION_PARAMETERS,DATA.DeltaT) ;
                VAR = KineticAndStrainEnergyCABLES(VAR,OPERFE,DATA,MATPRO) ;
            end
            break
        else
            %  Updating d_l by solving  Delta_d_L =  -K(DOFl,DOFl)\Res(DOFl)
            % 8. Assembly Stiffness Matrix (static part)
            
            K = GetStiffnessMatrixCABLES(OPERFE,DATA,VAR,FgradST,celastST) ;
            if DATA.ISDYNAMIC == 1
                % Dynamic part of the Jacobian
                K  =  KstifflDynamicPart(VAR,OPERFE,DATA.TIME_INTEGRATION_PARAMETERS,DATA.DeltaT,K) ;
            end
            
%             if ~isempty(OPERFE.HYDRO)
%                 % Contribution of hydrostatic forces
%                 % \FextPR{}{}   = - \NbST{T}  \wSTft{}
%                 % diag(kPRESSst)*Lbool
%                 % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/DynamicFloatingBodies.tex
%                 % ------------------------
%                 K = K - OPERFE.HYDRO.NbST_w'*(ConvertBlockDiag_general(kPRESSst,DATA.MESH.ndim,OPERFE.HYDRO.irowsNDIM,OPERFE.HYDRO.icolsNDIM)*OPERFE.HYDRO.Lbool) ;
%             end
            
            if DATA.SOLVER_IS_ITERATIVE == 1
                tolCONJG =[] ;
                niterCONJG =[] ;
                delta_dL=    pcg(K(DOFl,DOFl),-VAR.RESID(DOFl),[],[]) ;
            else
                delta_dL = -K(DOFl,DOFl)\VAR.RESID(DOFl);
            end
            VAR.DISP(DOFl) =  VAR.DISP(DOFl) + delta_dL ;
            normd = norm(delta_dL) ;
            
        end
        kiter = kiter + 1;
    end
end