function [DATA,CONVERGED]=NONLINEARdynamicCABLES(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO)

if nargin == 0
    load('tmp.mat')
end

% Initial acceleration
% -------------------
DATA = DefaultField(DATA,'InitialAccelerationFromResidual',0) ;
OPERFE.KinternalFORCES_given = []; 

if DATA.ISDYNAMIC  == 1 && DATA.InitialAccelerationFromResidual == 1
    VAR = IniAccelFromResidualGET(DATA,DISP_CONDITIONS,OPERFE,Fbody,Ftrac,VAR,MATPRO) ;
end

istep = 2;  icluster = 1;
nsteps = length(DATA.STEPS) ;
disp('*************************************')
disp(['LOOP over time steps (NSTEPS =',num2str(nsteps),')'])
disp('*************************************')
DOFr = DISP_CONDITIONS.DOFr ; DOFl = DISP_CONDITIONS.DOFl;
DATA = DefaultField(DATA,'SOLVER_IS_ITERATIVE',0) ;

while istep <=nsteps
    fprintf(1,'--------------------------------------------------\n')
    fprintf(1,['istep = %d      t_i = %4.5f \n'],[istep,DATA.STEPS(istep)])
    fprintf(1,'--------------------------------------------------\n')
    % 1.a) Initialization displacements
    DeltaT = DATA.STEPS(istep)-DATA.STEPS(istep-1) ;
    
    % For printing iterations --- We store as first step in this snapshot
    % matrices the converged at the previous time step
    iterk = 1;
    DATA = UpdateIterativeSnapshot(DATA,VAR,iterk) ;
    % ---------------------------------------------------------

    if DATA.ISDYNAMIC  == 1
        VAR = InertDampExternalForces(VAR,DATA.TIME_INTEGRATION_PARAMETERS,DeltaT,OPERFE) ;
        VAR.DISP(DOFl) = VAR.DISPtilde(DOFl); % Tilde variables ...
    end
    % Prescribed displacement at  time t_n+1
    if ~isempty(DISP_CONDITIONS.dR.U)
        dR = DISP_CONDITIONS.dR.U*DISP_CONDITIONS.dR.a(:,istep);
        VAR.DISP(DOFr) = dR;
    end
    % 1.b) Nodal  External forces (non-follower)
    VAR.FEXT = Fbody.U*Fbody.a(:,istep) + Ftrac.U*Ftrac.a(:,istep) ;
        
    % *********************************
    %%%% NEWTON-RAPHSON ALGORITHM
    % *****************************
    
    DATA.DeltaT = DeltaT ; DATA.istep = istep  ;
    [VAR,CONVERGED,DATA ]= NewtonRapshonStaticCABLES(DATA,OPERFE,VAR,MATPRO,DOFl) ;
        
    % Store snaphots    
    [SNAP,DATA,icluster] = StoreInfoSnapshots(istep,icluster,SNAP,VAR,DATA,CONVERGED) ;
    if CONVERGED == 0        
        disp('Convergence error ....')      
        break
    end    
    istep = istep + 1;
end