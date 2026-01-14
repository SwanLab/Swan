function [DATA,CONVERGED,QrotTIME,QrotINI]=NONLINEARstaticLARGE_COROT(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO)
% Adaptation of NONLINEARstaticLARGE.m to CO-rotational approach (EIFEM)
% JAHO, 29-oct-2024, Tuesday, Green's Aribau, Barcelona.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
% ------------------------------------------------------------

if nargin == 0
    load('tmp1.mat')
end
istep = 2;  icluster = 1;
nsteps = length(DATA.STEPS) ;
disp('*************************************')
disp(['LOOP over time steps (NSTEPS =',num2str(nsteps),')'])
disp('*************************************')
DOFr = DISP_CONDITIONS.DOFr ; DOFl = DISP_CONDITIONS.DOFl;
if isfield(DISP_CONDITIONS,'DOFm')
    DOFm = DISP_CONDITIONS.DOFm ;
    OPERFE.A = DISP_CONDITIONS.A;
    OPERFE.G = DISP_CONDITIONS.G;
    OPERFE.DOFr =  DOFr ;
    OPERFE.DOFm = DOFm ;
else
    DOFm = [] ;
    OPERFE.A = [];
end
DATA = DefaultField(DATA,'SOLVER_IS_ITERATIVE',0) ;
DATA = DefaultField(DATA,'SMALL_STRAIN_KINEMATICS',0) ;
DATA.NEWTON_RAPHSON = DefaultField(DATA.NEWTON_RAPHSON,'ENFORCE_BOTH_TOLERANCES_RESIDUAL_AND_DISPLACEMENTS',0) ;
OPERFE = DefaultField(OPERFE,'KinternalFORCES_given',[])  ;  % In case STIFFNESS MATRIX is given by the user
OPERFE = DefaultField(OPERFE,'DOFm',[])  ;  % 
DATA = DefaultField(DATA,'CECM_ONLY_FOR_NONLINEAR_STRESSES',0) ; % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
DATA.MESH = DefaultField(DATA.MESH,'ndimFINE',DATA.MESH.ndim) ;
DATA.MESH.ndim = DATA.MESH.ndimFINE;  

% ROTATION MATRICES, Initialization
% ----------------------------------------
% \DiagC{\QrotALL} = \DiagC{\QrotINIall}
Qrot = OPERFE.QrotINI; 
QrotINI = OPERFE.QrotINI ; 

%
TOL_CONVERGENCE_angle_rads = 1e-3;  
DATA = DefaultField(DATA,'TOL_CONVERGENCE_angle_rads',TOL_CONVERGENCE_angle_rads) ;
NITER_ROTATIONS_UPDATE= 10 ; 
DATA = DefaultField(DATA,'NITER_ROTATIONS_UPDATE',NITER_ROTATIONS_UPDATE) ;

 QrotTIME = cell(1,nsteps) ; 
 
% profile on 
while istep <=nsteps
    fprintf(1,'--------------------------------------------------\n')
    fprintf(1,['istep = %d      t_i = %4.5f \n'],[istep,DATA.STEPS(istep)])
    fprintf(1,'--------------------------------------------------\n')
    % 1.a) Initialization displacements
    
    % Prescribed macroscopic deformation at  time t_n+1
    dR = DISP_CONDITIONS.dR.U*DISP_CONDITIONS.dR.a(:,istep);
    OPERFE.uBAR = dR;
    if ~isempty(DOFm)
        % See explanation in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/10_PeriodicQUADL.mlx
        % \d^{n+1,k}_L =   \d^{n}_L
        % \d^{n+1,k}_S  = \A  \d^{n+1,k}_M + \uBAR(t)
        % Affine boundary conditions
        dR = dR + DISP_CONDITIONS.G*VAR.DISP(DOFm);
        %VAR.DISP(DOFr) = dR;
    end
    VAR.DISP(DOFr) = dR;
    % 1.b) External forces
    if ~isempty(Fbody.U)
        VAR.FEXT = Fbody.U*Fbody.a(:,istep)  ; % change, 29-Jan-2023
    end
    
    if ~isempty(Ftrac.U)
        VAR.FEXT = VAR.FEXT + Ftrac.U*Ftrac.a(:,istep) ; % Decomposition Space time (via SVD)
    end
    
    
    % *********************************
    %%%% NEWTON-RAPHSON ALGORITHM + UPDATE ROTATIONS 
    % *****************************+
    iterROTATIONS = 1;
    norm_angle_rads = 1e20 ; 
    while iterROTATIONS <= DATA.NITER_ROTATIONS_UPDATE && norm_angle_rads > DATA.TOL_CONVERGENCE_angle_rads
        [BstQ,KcLINq,dCqLOC,D_QrotALL] = UpdateMatricesRotation(OPERFE,Qrot) ;
        [VAR,CONVERGED ]= NewtonRapshonStaticLarge_COROT(DATA,OPERFE,VAR,MATPRO,DOFl,dCqLOC,KcLINq,BstQ,D_QrotALL) ;
        % COMPUTE ROTATION ANGLES and update 
        [Qrot,norm_angle_rads] = UpdateRotationsEIFEM(OPERFE,Qrot,VAR,DATA,iterROTATIONS) ;
        iterROTATIONS = iterROTATIONS +1 ; 
    end
    
    QrotTIME{istep} = Qrot ; 
    % Store snaphots
    
    [SNAP,DATA,icluster] = StoreInfoSnapshots(istep,icluster,SNAP,VAR,DATA,CONVERGED) ;
    if CONVERGED == 0 || norm_angle_rads > DATA.TOL_CONVERGENCE_angle_rads
        
        if CONVERGED == 0
        disp('Convergence error (lack of equilibrium, residual forces) ....')
        % load gong.mat;
        %  soundsc(y)
        pause
        break
        else
          disp('Warning: Unable to find intermediate configuration (co-rotational)')  
            
        end
        
        
        
    end
    
    istep = istep + 1;
end


%profile report