function [DATA,CONVERGED]=NONLINEARstaticCABLES(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO)

if nargin == 0
    load('tmp.mat')
end

istep = 2;  icluster = 1;
nsteps = length(DATA.STEPS) ;
disp('*************************************')
disp(['LOOP over time steps (NSTEPS =',num2str(nsteps),')'])
disp('*************************************')
DOFr = DISP_CONDITIONS.DOFr ; DOFl = DISP_CONDITIONS.DOFl;
DATA = DefaultField(DATA,'SOLVER_IS_ITERATIVE',0) ;
DATA.NEWTON_RAPHSON = DefaultField(DATA.NEWTON_RAPHSON,'ENFORCE_BOTH_TOLERANCES_RESIDUAL_AND_DISPLACEMENTS',0) ;
%OPERFE = DefaultField(OPERFE,'KinternalFORCES_given',[])  ;  % In case STIFFNESS MATRIX is given by the user 


while istep <=nsteps
    fprintf(1,'--------------------------------------------------\n')
    fprintf(1,['istep = %d      t_i = %4.5f \n'],[istep,DATA.STEPS(istep)])
    fprintf(1,'--------------------------------------------------\n')
    % 1.a) Initialization displacements
    
    % Prescribed macroscopic deformation at  time t_n+1
    dR = DISP_CONDITIONS.dR.U*DISP_CONDITIONS.dR.a(:,istep);
    VAR.DISP(DOFr) = dR;
    % 1.b) External forces
    VAR.FEXT = Fbody.U*Fbody.a(:,istep) + Ftrac.U*Ftrac.a(:,istep) ;
    
    % *********************************
    %%%% NEWTON-RAPHSON ALGORITHM
    % *****************************+
      if istep ==26
          disp('borrar estooo')
      end
    
    [VAR,CONVERGED,DATA]= NewtonRapshonStaticCABLES(DATA,OPERFE,VAR,MATPRO,DOFl) ;
    
    % Store snaphots
    
    [SNAP,DATA,icluster] = StoreInfoSnapshots(istep,icluster,SNAP,VAR,DATA,CONVERGED) ;
    if CONVERGED == 0
        
        disp('Convergence error ....')
       % load gong.mat;
      %  soundsc(y)
        pause
        break
    end
    
    istep = istep + 1;
end