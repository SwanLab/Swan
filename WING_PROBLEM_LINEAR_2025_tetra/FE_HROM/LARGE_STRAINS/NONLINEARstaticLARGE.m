function [DATA,CONVERGED]=NONLINEARstaticLARGE(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO)
% ---------------------------------------------------------------------------------------------------
% FUNCTION: NONLINEARstaticLARGE
% ---------------------------------------------------------------------------------------------------
% PURPOSE:
%   This function performs a nonlinear static analysis over multiple time steps using the
%   Newton–Raphson method. It handles general displacement-controlled boundary conditions,
%   body and traction forces, and updates the solution history (snapshots) for postprocessing
%   or model reduction purposes.
%
% USAGE:
%   [DATA, CONVERGED] = NONLINEARstaticLARGE(DATA, DISP_CONDITIONS, VAR, OPERFE, ...
%                                            SNAP, Fbody, Ftrac, MATPRO)
%
% INPUTS:
%   - DATA            : Structure with simulation parameters, time discretization, and solver settings.
%   - DISP_CONDITIONS : Structure with prescribed displacements, boundary condition mappings
%                       (DOFl, DOFr, possibly DOFm for affine conditions), and time-dependent control vectors.
%   - VAR             : Structure containing solution fields (displacements, forces, internal variables).
%   - OPERFE          : Structure with finite element operators (e.g., A, G, mappings, preassembled stiffness).
%   - SNAP            : Structure for storing snapshots of key variables across time steps.
%   - Fbody           : Structure containing body force loadings (matrices `U` and coefficients `a`).
%   - Ftrac           : Structure with surface traction loadings (matrices `U` and coefficients `a`).
%   - MATPRO          : Material property structure used in the constitutive model.
%
% OUTPUTS:
%   - DATA        : Updated simulation data structure (time, convergence, snapshot control, etc.).
%   - CONVERGED   : Flag (1 if converged at current step, 0 if not).
%
% FUNCTIONALITY:
%   - Initializes the displacements at each time step based on Dirichlet boundary conditions.
%   - Applies affine boundary conditions using projection matrices (A, G, uBAR) if required.
%   - Assembles external forces due to body and surface loads (including SVD-decomposed representations).
%   - Calls the Newton–Raphson nonlinear solver at each time step (`NewtonRapshonStaticLarge`).
%   - Stores snapshots of displacements, internal variables, and other fields for further analysis.
%   - Detects and stops upon convergence failure.
%
% FEATURES:
%   - Compatible with periodic/multiscale setups using affine constraints.
%   - Modular treatment of external loadings and boundary conditions.
%   - Uses optional solver configurations like exact/approximate Jacobians, user-defined stiffness,
%     or convergence override in the presence of negative Jacobians.
%
% NOTES:
%   - This function is time-stepping oriented but tailored to quasi-static (or dynamic via inner flags).
%   - Snapshots can be used later for reduced-order modeling or postprocessing in GiD.
%
% AUTHOR:
%   Joaquín A. Hernández, UPC – CIMNE
%   Date: 29-Jan-2023 (last updated for space–time decomposition features)
%   Comments by ChatGPT-4, on 12-May-2025
% DEPENDENCIES:
%   - NewtonRapshonStaticLarge
%   - StoreInfoSnapshots
%   - DefaultField (utility to handle optional fields)
%
% ---------------------------------------------------------------------------------------------------

if nargin == 0
    load('tmp.mat')
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
    
    
    % 1.b) External forces
    if ~isempty(Fbody.U)
        VAR.FEXT = Fbody.U*Fbody.a(:,istep)  ; % change, 29-Jan-2023
    end
    
    if ~isempty(Ftrac.U)
        VAR.FEXT = VAR.FEXT + Ftrac.U*Ftrac.a(:,istep) ; % Decomposition Space time (via SVD)
    end
    
    
    % *********************************
    %%%% NEWTON-RAPHSON ALGORITHM
    % *****************************+
    
    [VAR,CONVERGED ]= NewtonRapshonStaticLarge(DATA,OPERFE,VAR,MATPRO,DOFl) ;
    
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