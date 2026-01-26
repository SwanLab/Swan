function [DATA,CONVERGED,MACROVAR]=NONLINEARhomogLARGE(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,MATPRO)

if nargin == 0
    load('tmp1.mat')
end

istep = 2;  icluster = 1;
nsteps = length(DATA.STEPS) ;
disp('*************************************')
disp(['LOOP over time steps (NSTEPS =',num2str(nsteps),')'])
disp('*************************************')
DOFr = DISP_CONDITIONS.DOFr ; DOFl = DISP_CONDITIONS.DOFl; DOFm = DISP_CONDITIONS.DOFm;
DATA = DefaultField(DATA,'SOLVER_IS_ITERATIVE',0) ; 

% MACROVAR.PK1STRESS = zeros(DATA.MESH.ndim^2,nsteps) ; 
% MACROVAR.PK2STRESS = zeros(DATA.MESH.nstrain,nsteps) ; 

while istep <=nsteps
    disp('--------------------------------------------------')
    disp(['istep = ',num2str(istep), '  t_i =',num2str(DATA.STEPS(istep))])
    disp('---------------------------------------------------')
    % 1.a) Initialization displacements
    
    % Prescribed macroscopic deformation gradeint at  time t_n+1
    GRADuMACROst = DISP_CONDITIONS.MACROVAR.GRADuMACROst.U*DISP_CONDITIONS.MACROVAR.GRADuMACROst.a(:,istep) ; 
    kiter = 1;
    CONVERGED = 0 ;
    normd = 1e20 ;
    while  kiter <=DATA.NEWTON_RAPHSON.NMAXiter
        % 2. Deformation gradient at all Gauss points
        FgradST = GRADuMACROst + OPERFE.BstA*VAR.FLUCT(DOFl)  + OPERFE.IDENTITY_F  ;
        % 3. Green-Lagrante strains at all Gauss points
        VAR.GLSTRAINS = StrainGreenLagrange(FgradST,DATA.MESH.ndim) ;
        % 4. 2nd Piola-Kirchhoff stresses at all Gauss Points
        [VAR.PK2STRESS,celastST ]= PK2stress_Constitutive_Model(VAR.GLSTRAINS,MATPRO,DATA,FgradST) ;
        % 5. 1st Piola-Kirchhoff stresses at all Gauss Points
        PoneST = PK1stress(VAR.PK2STRESS,FgradST,DATA.MESH.ndim) ;
        % 6.1. Internal forces
        VAR.RESID  = InternalForcesHOMOG(OPERFE,PoneST,DATA) ; 
        % 7. Check convergence
        [criter_f, norm_res]=  CheckConvergenceHOMOG(VAR.RESID,kiter,normd) ;
        
        if criter_f <= DATA.NEWTON_RAPHSON.TOL_FORCES_REL ||  normd <= DATA.NEWTON_RAPHSON.TOL_DISPLAC_ABS
            % Converged
            CONVERGED = 1;
            % Compute Von Mises Stress
            [ VAR.VONMISES_CAUCHY_STRESS,VAR.CAUCHY_STRESS ]= ...
                VonMisesCauchyStresses(PoneST,FgradST,DATA.MESH.ndim,[],DATA) ;
            if ~isempty(DOFr)
            VAR.FLUCT(DOFr) = DISP_CONDITIONS.G*VAR.FLUCT(DOFm) ; 
            end
            
            % UPDATE TOTAL DISPLACEMENTS 
          %  VAR.DISP =  DISP_CONDITIONS.MACROVAR.DISP.U* DISP_CONDITIONS.MACROVAR.DISP.a(:,istep) + VAR.FLUCT ; 
            VAR.DISP =  VAR.FLUCT ; % The macro-displacement is addedd when printing in GID 
            
            % Compute homogeneized stresses 
            for istrain = 1:DATA.MESH.ndim^2
                IND = istrain:DATA.MESH.ndim^2:size(PoneST,1) ; 
                VAR.PK1STRESS_homog(istrain)= OPERFE.WEIGHTShomog'*PoneST(IND);
            end
            
%             for istrain = 1:DATA.MESH.nstrain
%                 IND = istrain:DATA.MESH.nstrain:size(VAR.PK2STRESS,1) ; 
%                 MACROVAR.PK2STRESS(istrain,istep) = OPERFE.WEIGHTShomog'*VAR.PK2STRESS(IND);
%             end
            
            
            break
        else
            %  Updating d_l by solving  Delta_d_L =  -K(DOFl,DOFl)\Res(DOFl)
            % 8. Assembly Stiffness Matrix
            K = KstiffLargeStrainsHOMOG(OPERFE,VAR.PK2STRESS,FgradST,DATA.MESH.ndim,celastST) ;
            if DATA.SOLVER_IS_ITERATIVE == 1
                tolCONJG =[] ;
                niterCONJG =[] ;
                error('Revised this part')
                delta_dL=    pcg(K(DOFl,DOFl),-VAR.RESID(DOFl),[],[]) ;
            else
                delta_dL = -K\VAR.RESID;
            end
            VAR.FLUCT(DOFl)=  VAR.FLUCT(DOFl) + delta_dL ;
            normd = norm(delta_dL) ;
        end
        kiter = kiter + 1;
    end
    
    
    [SNAP,DATA,icluster] = StoreInfoSnapshots(istep,icluster,SNAP,VAR,DATA,CONVERGED) ;
    if CONVERGED == 0
        
        disp('Convergence error ....')
        break
    end
    
    istep = istep + 1;
end