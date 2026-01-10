function [DATA,CONVERGED]=NONLINEARstaticLARGE_manifoldFST(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO)
%--------------------------------------------------------------------------
% NONLINEARstaticLARGE_manifoldFST
%
% PURPOSE:
%   Faster, updated variant of `NONLINEARstaticLARGE_manifold` that performs
%   a sequence of nonlinear STATIC solves over prescribed time steps using a
%   manifold-based ROM. This FST version is compatible with the FAST pipeline
%   (no anonymous-function evaluators for τ(q)); when available, it leverages a
%   precompiled evaluator structure for the manifold-ECM reconstruction of
%   internal work densities (master→slave) during Newton iterations.
%
% WHAT THIS ROUTINE DOES (per time step):
%   1) Enforces Dirichlet boundary conditions (standard or affine/periodic via A,G).
%   2) Assembles external loads from space–time separated forms (Fbody, Ftrac).
%   3) Solves the static equilibrium using a Newton–Raphson scheme on the
%      tangent space of the manifold:
%         • If no manifold-ECM evaluator is present:
%               `NewtonRapshonStaticLarge_manifoldFST`
%         • If manifold-ECM evaluator exists (OPERFE.DATA_regress_eta_der):
%               `NewtonRapshonStaticLarge_manifoldECMnonFST`
%      (the latter reconstructs slave internal-work contributions from the
%       master point, reducing cost while preserving accuracy).
%   4) Stores snapshots (displacements, forces, etc.) via `StoreInfoSnapshots`.
%
% KEY DIFFERENCES VS. ORIGINAL (`NONLINEARstaticLARGE_manifold`):
%   • Removes reliance on τ(q), τ′(q), τ″(q) anonymous-function handles.
%   • When available, uses `OPERFE.DATA_regress_eta_der` (produced offline by
%     `DiscreteECM_adaptWEIGHTSfst` → `BsplinesLeastSquares_fastECM`) to compute
%     modal amplitudes / internal-work reconstructions efficiently.
%   • Selects the appropriate Newton solver automatically based on the presence
%     of manifold-ECM data.
%
% INPUTS:
%   - DATA            : Main simulation struct (steps, tolerances, flags, mesh dims).
%   - DISP_CONDITIONS : Partition of DOFs (DOFl, DOFr, optional DOFm) and BC data
%                       (dR.U, dR.a, optional A,G for affine/periodic constraints).
%   - VAR             : State container (displacements, residuals, etc.) at entry time.
%   - OPERFE          : Reduced FE operators and auxiliary maps:
%                         • Bst, M, (optional) KinternalFORCES_given
%                         • uBAR, DOFl/DOFr/DOFm
%                         • (optional) DATA_regress_eta_der for manifold-ECM
%   - SNAP            : Snapshot accumulator/metadata (cluster info, files).
%   - Fbody, Ftrac    : Space–time separated external forces (U, a(t)).
%   - MATPRO          : Material model parameters / constitutive data.
%
% OUTPUTS:
%   - DATA       : Updated simulation struct (convergence flags, snapshot info).
%   - CONVERGED  : Logical flag (1 if all steps converged, 0 otherwise).
%
% THEORY / CONTEXT:
%   - Manifold ROM with tangent projection (τ(q)) and ECM hyperreduction, as in
%     MLEARNstruct_1.pdf:
%       §9.2  – Manifold mapping and tangent-space Newton iterations.
%       §9.3  – ECM and master/slave reconstruction of internal-work densities.
%
% DEPENDENCIES:
%   - NewtonRapshonStaticLarge_manifoldFST
%   - NewtonRapshonStaticLarge_manifoldECMnonFST
%   - StoreInfoSnapshots
%   - DefaultField
%
% HISTORY (CREATION / UPDATES):
%   • Created  : 22-Jul-2025 (as `NONLINEARstaticLARGE_manifold`)
%   • Updated  : 11-Aug-2025, Molinos Marfagones (Cartagena)
%                - Introduced FAST path: switched to precompiled manifold-ECM
%                  evaluator (`OPERFE.DATA_regress_eta_der`) and added automatic
%                  solver selection (standard vs ECM-nonlinear).
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC/CIMNE)
%   Comments updated by ChatGPT, 12-Aug-2025
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
DATA.MESH = DefaultField(DATA.MESH,'ndimFINE',DATA.MESH.ndim) ;
DATA.MESH.ndim = DATA.MESH.ndimFINE;

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
    
    if  isempty(OPERFE.DATA_regress_eta_der)
        % Standard ECM 
        if isempty(OPERFE.wSTs_cluster)
          [VAR,CONVERGED ]= NewtonRapshonStaticLarge_manifoldFST(DATA,OPERFE,VAR,MATPRO,DOFl) ;
        else
            % VAriable weights 
           [VAR,CONVERGED ]= NewtonRapshonStaticLarge_manifoldFSTw(DATA,OPERFE,VAR,MATPRO,DOFl) ;  
        end
    else
        % New version (21-July-2025), see
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/03_ECMmanifold.mlx
        % It leverages the nonlinear manifold in which internal work
        % density lives 
        [VAR,CONVERGED ]= NewtonRapshonStaticLarge_manifoldECMnonFST(DATA,OPERFE,VAR,MATPRO,DOFl) ;
    end
    
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