function [DATA,CONVERGED]=NONLINEARstaticLARGE_bubfe(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/08_IndividualTRAINbub.mlx%
% Self-equilibrium problem BUBBLE APPROACH
% JAHO, 30-OCT-2023, BARCELONA, BALMES 185
% --------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end
istep = 2;  icluster = 1;
nsteps = length(DATA.STEPS) ;
disp('*************************************')
disp(['LOOP over time steps (NSTEPS =',num2str(nsteps),')'])
disp('*************************************')
%DOFr = DISP_CONDITIONS.DOFr ; DOFl = DISP_CONDITIONS.DOFl;
DATA = DefaultField(DATA,'SOLVER_IS_ITERATIVE',0) ;
DATA = DefaultField(DATA,'SMALL_STRAIN_KINEMATICS',0) ;
DATA.NEWTON_RAPHSON = DefaultField(DATA.NEWTON_RAPHSON,'ENFORCE_BOTH_TOLERANCES_RESIDUAL_AND_DISPLACEMENTS',0) ;
DATA.FirstIterationAtWhichYieldingIsActivated = 0 ; % this is a local amendment for solving the issue appearing in 
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/FE_CODE/NONLINEAR_Dyn/SmallStrainJ2Plasticity_LOCAL.m
%OPERFE = DefaultField(OPERFE,'KinternalFORCES_given',[])  ;  % In case STIFFNESS MATRIX is given by the user 
OPERFE.BstA = OPERFE.Bst*DISP_CONDITIONS.A; 


while istep <=nsteps
    fprintf(1,'--------------------------------------------------\n')
    fprintf(1,['istep = %d      t_i = %4.5f \n'],[istep,DATA.STEPS(istep)])
    fprintf(1,'--------------------------------------------------\n')
    % 1.a) Initialization bubb
    
    % Prescribed coarse-scale displacements at step istep 
    VAR.ELAST_DISP = DISP_CONDITIONS.PhiDEF*DISP_CONDITIONS.COARSE_SCALE_HISTORY(:,istep)  ;
    VAR.DISP = VAR.ELAST_DISP + VAR.BUB_DISP ; 
    
    
%     dR = DISP_CONDITIONS.dR.U*DISP_CONDITIONS.dR.a(:,istep);
%     VAR.DISP(DOFr) = dR;
%     % 1.b) External forces
%     VAR.FEXT = Fbody.U*Fbody.a(:,istep)  ; % change, 29-Jan-2023
%     if ~isempty(Ftrac.U)
%         VAR.FEXT = VAR.FEXT + Ftrac.U*Ftrac.a(:,istep) ;
%     end
     
    
    % *********************************
    %%%% NEWTON-RAPHSON ALGORITHM
    % *****************************+
    
    [VAR,CONVERGED ]= NewtonRapshonStaticLarge_bub(DATA,OPERFE,VAR,MATPRO,DISP_CONDITIONS) ;
    
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