function [DATA,CONVERGED]=NONLINEARstaticLARGE_manifoldMD(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO)
%--------------------------------------------------------------------------
% NONLINEARstaticLARGE_manifoldMD
% --------------------------------
% DESCRIPTION
%   Static manifold-ROM driver for geometrically nonlinear problems over a
%   prescribed sequence of “time” steps (used here as load increments). This
%   routine generalizes NONLINEARstaticLARGE_manifoldFST to support more than
%   one generalized coordinate, enforces affine (A,G) or standard Dirichlet
%   boundary conditions, assembles space–time separated loads, and advances
%   the solution with Newton–Raphson on the tangent manifold. It also handles
%   adaptive solver safeguards (geometric-term disabling, optional line search)
%   and snapshot storage for post-processing.
%
% WHAT IT DOES (per increment i)
%   1) Boundary conditions:
%      - Build prescribed displacement vector dR = dR.U * dR.a(:,i).
%      - If periodic/affine BCs exist (G, A, DOFm), shift dR ← dR + G * u_m
%        (with u_m from the current state VAR.DISP(DOFm)), assign DOFr values.
%   2) External actions:
%      - Assemble space–time separated loads:
%           FEXT = Fbody.U * Fbody.a(:,i)  +  Ftrac.U * Ftrac.a(:,i).
%   3) Nonlinear solve (Newton–Raphson on the manifold):
%      - If no manifold-ECM evaluator is available (OPERFE.DATA_regress_eta_der empty):
%          • If no ECM clusters (OPERFE.wSTs_cluster empty):
%              Try with geometric stiffness term (DATA.USE_GEOMETRIC_term_K = 1);
%              on failure, retry with geometric term disabled (=0).
%          • Else (clusters present):
%              Use either NewtonRaphsonStrategy_disablingTERMS or, if enabled,
%              a line-search variant (NewtonRapshonStatic_MHROM_LSearch_GROK).
%      - If manifold-ECM evaluator exists:
%          Path currently abandoned (October 2025) — call raises an error by design.
%   4) Snapshots:
%      - StoreInfoSnapshots saves states/metadata; convergence flag is recorded.
%
% KEY FEATURES / SAFEGUARDS
%   • Multi-DOF manifold support (more than one generalized coordinate).
%   • Affine/periodic BC handling via matrices (A,G) and DOFm partition.
%   • FEXT assembled from space–time separated factors (U, a(t)) for Fbody/Ftrac.
%   • Automatic fallback: if the geometric term hinders convergence, it is
%     disabled and the Newton solve is retried.
%   • Optional line search (DATA.Newton_Raphson_WithLineSearch = true) routed to
%     NewtonRapshonStatic_MHROM_LSearch_GROK.
%   • Snapshot accumulation via StoreInfoSnapshots with cluster bookkeeping.
%
% INPUTS (summary)
%   DATA            : Global run config (STEPS, NEWTON_RAPHSON tolerances,
%                     flags like SMALL_STRAIN_KINEMATICS, PROCEED_WITH_NEGATIVE_JACOBIANS,
%                     DisableGeometricAndWeightTermWhenNoConvergence, etc.).
%   DISP_CONDITIONS : DOF partitions (DOFl, DOFr, optional DOFm) and BC data
%                     (dR.U, dR.a, optional A,G for affine/periodic).
%   VAR             : State at entry (e.g., VAR.DISP, residuals, etc.).
%   OPERFE          : Reduced operators/maps (Bst, M, optional KinternalFORCES_given,
%                     optional DATA_regress_eta_der for manifold-ECM, uBAR holder).
%   SNAP            : Snapshot accumulator (cluster indices, file names).
%   Fbody, Ftrac    : Space–time separated loads with fields .U and .a.
%   MATPRO          : Constitutive parameters / model data.
%
% OUTPUTS
%   DATA      : Updated configuration (convergence info, snapshot updates).
%   CONVERGED : Logical (true if all increments solved; false if an increment failed).
%
% DEPENDENCIES (called here)
%   - NewtonRapshonStaticLarge_manifoldMD
%   - NewtonRaphsonStrategy_disablingTERMS
%   - NewtonRapshonStatic_MHROM_LSearch_GROK
%   - StoreInfoSnapshots
%   - DefaultField
%   (abandoned path) NewtonRapshonStaticLarge_manifoldECMnonMD  [raises error]
%
% HISTORY (keep existing references)
%   - Generalizes NONLINEARstaticLARGE_manifoldFST.
%   - Date modification: 13th August 2025, Wednesday — Molinos Marfagones, Cartagena.
%
% AUTHOR
%   Joaquín A. Hernández Ortega (UPC/CIMNE)
%
% COMMENTS UPDATED BY: ChatGPT-5
%   Date of update: 07-NOV-2025, Barcelona, Spain
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

% ---------------------------------------------------------------------------------------------------

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
DATA = DefaultField(DATA,'PROCEED_WITH_NEGATIVE_JACOBIANS',1) ;
DATA = DefaultField(DATA,'DisableGeometricAndWeightTermWhenNoConvergence',true) ;
DATA = DefaultField(DATA,'Newton_Raphson_WithLineSearch',false) ;
DATA = DefaultField(DATA,'FreezeBasis',true) ;


DATA.MESH = DefaultField(DATA.MESH,'ndimFINE',DATA.MESH.ndim) ;
DATA.MESH.ndim = DATA.MESH.ndimFINE;

while istep <=nsteps
    fprintf(1,'--------------------------------------------------\n')
    fprintf(1,['istep = %d      t_i = %4.5f \n'],[istep,DATA.STEPS(istep)])
    fprintf(1,'--------------------------------------------------\n')
    % 1.a) Initialization displacements
    
    %
%      if istep == 13
%         warning('BORRAR ESTO')
%     end
%     
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
    OPERFE.DOFr = DOFr;
    OPERFE.DOFl = DOFl;
    % 1.b) External forces
    if ~isempty(Fbody.U)
        VAR.FEXT_extended = Fbody.U*Fbody.a(:,istep)  ; % change, 29-Jan-2023
    end
    
    if ~isempty(Ftrac.U)
        VAR.FEXT_extended = VAR.FEXT_extended + Ftrac.U*Ftrac.a(:,istep) ; % Decomposition Space time (via SVD)
    end
    
    
    % *********************************
    %%%% NEWTON-RAPHSON ALGORITHM
    % *****************************+
    
    if  isempty(OPERFE.DATA_regress_eta_der) % || isempty(OPERFE.wSTs_cluster)
        % Standard ECM
        if isempty(OPERFE.wSTs_cluster)
            
          %  DATA = DefaultField(DATA,'FreezeLeftBasis',true) ;
%            if DATA.FreezeLeftBasis == 1 
%                DATA.FreezeLeftBasis = 1;   
%                 DATA.USE_GEOMETRIC_term_K = 0; 
%                 [VAR,CONVERGED ]= NewtonRapshonStaticLarge_manifoldMD(DATA,OPERFE,VAR,MATPRO,DOFl) ;
%            else
            DATA.USE_GEOMETRIC_term_K = 1;            
            [VAR,CONVERGED ]= NewtonRapshonStaticLarge_manifoldMD(DATA,OPERFE,VAR,MATPRO,DOFl) ;
            
 
            
            
            if CONVERGED ==0
                disp('Disabling geometric term')
                DATA.USE_GEOMETRIC_term_K = 0;                
                [VAR,CONVERGED ]= NewtonRapshonStaticLarge_manifoldMD(DATA,OPERFE,VAR,MATPRO,DOFl) ;
            end
       %    end
            
            
            
        else
            
            if DATA.Newton_Raphson_WithLineSearch
               %  [VAR,CONVERGED ]= NewtonRapshonStatic_MHROM_LSearch(DATA,OPERFE,VAR,MATPRO,DOFl) ;
                 [VAR,CONVERGED ]= NewtonRapshonStatic_MHROM_LSearch_GROK(DATA,OPERFE,VAR,MATPRO,DOFl) ;
      %    [VAR,CONVERGED ]= NewtonRapshonStatic_MHROM_LSearch_GROK_disab(DATA,OPERFE,VAR,MATPRO,DOFl) ;
    
                 
          %       [VAR,CONVERGED ]= NewtonRapshonStatic_MHROM_LSearch_GROK_v1(DATA,OPERFE,VAR,MATPRO,DOFl) ;

            %  [VAR,CONVERGED ]= NewtonRapshonStatic_MHROM_LSearch_CGPTv1(DATA,OPERFE,VAR,MATPRO,DOFl) ;
      % [VAR,CONVERGED ]= NewtonRapshonStatic_MHROM_LSearch_CGPTv2(DATA,OPERFE,VAR,MATPRO,DOFl) ;


            else
                [VAR,CONVERGED ]= NewtonRaphsonStrategy_disablingTERMS(DATA,OPERFE,VAR,MATPRO,DOFl) ;
            end
            
        end
        
    else
        % New version (21-July-2025), see
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/03_ECMmanifold.mlx
        % It leverages the nonlinear manifold in which internal work
        % density lives
        % Update for 2 latent variables in
        %  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/08_PLAST_adapECx
        error('Option abandoned...October 2025')
        [VAR,CONVERGED ]= NewtonRapshonStaticLarge_manifoldECMnonMD(DATA,OPERFE,VAR,MATPRO,DOFl) ;
        
    end
    
    % Store snaphots
    
    [SNAP,DATA,icluster] = StoreInfoSnapshots(istep,icluster,SNAP,VAR,DATA,CONVERGED) ;
    if CONVERGED == 0
        
        disp('Convergence error ....')
        % load gong.mat;
        %  soundsc(y)
        disp(['qPLAST = ',num2str(VAR.DISP(:)')])
        pause
        break
    end
    
    istep = istep + 1;
end