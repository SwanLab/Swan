function [DATA,CONVERGED]=NONLINEARstaticLARGE_manifold(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO)
%--------------------------------------------------------------------------
% function [DATA, CONVERGED] = ...
%     NONLINEARstaticLARGE_manifold(DATA, DISP_CONDITIONS, VAR, ...
%                                   OPERFE, SNAP, Fbody, Ftrac, MATPRO)
%
% PURPOSE:
%   Executes a nonlinear **static simulation** across a sequence of time steps
%   using a reduced-order model (ROM) in which the displacement field is 
%   parametrized over a nonlinear manifold τ(q). 
%
%   For each time step:
%     - Boundary conditions (affine or periodic) are enforced.
%     - External forces (e.g., body or traction) are applied via SVD-reduced bases.
%     - Newton–Raphson iterations are performed in the tangent space of the manifold.
%     - Snapshots of displacements, forces, and internal variables are stored.
%
%   Two Newton solvers are supported:
%     1. `NewtonRapshonStaticLarge_manifold` – standard projection with τ(q)
%     2. `NewtonRapshonStaticLarge_manifoldECMnon` – enhanced version that
%        accounts for hyperreduction via nonlinear internal force extrapolation (η_NON).
%
% INPUTS:
% --------
%   DATA : struct
%       Global simulation parameters, tolerances, flags, and time integration.
%
%   DISP_CONDITIONS : struct
%       Contains DOF partitioning (DOFl, DOFr, DOFm) and boundary conditions:
%       - uBAR: prescribed affine macro deformation
%       - A, G: mapping matrices for enforcing affine or periodic BCs
%
%   VAR : struct
%       Current solution fields (e.g., displacement, stress, residual).
%
%   OPERFE : struct
%       Reduced operators and manifold-related maps (τ, ∂τ/∂q, η_NON, etc.).
%
%   SNAP : struct
%       Container for snapshot storage across time steps.
%
%   Fbody, Ftrac : struct
%       Space-time decomposed external forces via reduced bases (e.g., from SVD).
%
%   MATPRO : struct
%       Material model definition and constitutive parameters.
%
% OUTPUTS:
% ---------
%   DATA : struct
%       Updated simulation structure, now including stored solutions and convergence info.
%
%   CONVERGED : logical
%       Flag (1 if converged at every time step, 0 otherwise).
%
% FEATURES:
% ---------
%   - Supports manifold-based hyperreduction using ECM and nonlinear extrapolation of
%     internal force densities (η_NON).
%
%   - Dynamically chooses the appropriate Newton solver depending on the presence
%     of `OPERFE.etaNON`.
%
%   - Handles affine or periodic boundary conditions using transformation matrices A, G.
%
%   - Collects snapshots at each time step for reconstruction, analysis, or machine learning.
%
% REFERENCES:
%   - Hernández Ortega, J.A., *Machine Learning Techniques in Structural Analysis*
%     Section 9.2 (pp. 71–77): Nonlinear manifold-based model reduction.
%
% DEPENDENCIES:
%   - NewtonRapshonStaticLarge_manifold
%   - NewtonRapshonStaticLarge_manifoldECMnon
%   - StoreInfoSnapshots
%   - DefaultField
%
% AUTHOR:
%   Joaquín A. Hernández Ortega, UPC – CIMNE
%   Last revised: 22 July 2025, Honest Greens – Pedralbes, Barcelona
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
    if  isempty(OPERFE.etaNON)
    [VAR,CONVERGED ]= NewtonRapshonStaticLarge_manifold(DATA,OPERFE,VAR,MATPRO,DOFl) ;
    else
        % New version (21-July-2025), see
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/03_ECMmanifold.mlx
        % It leverages the nonlinear manifold in which internal work
        % density lives 
        [VAR,CONVERGED ]= NewtonRapshonStaticLarge_manifoldECMnon(DATA,OPERFE,VAR,MATPRO,DOFl) ;
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