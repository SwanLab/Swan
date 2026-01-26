function [DATA,CONVERGED,QrotTIME]=NONLINEARstaticLARGE_COROT_LRss(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO)
% ---------------------------------------------------------------------------------------------------
% FUNCTION: NONLINEARstaticLARGE_COROT_LRss
% ---------------------------------------------------------------------------------------------------
% PURPOSE:
%   This function performs nonlinear static analysis for mechanical systems under **small strains
%   and large rotations** (or vice versa), using a **corotational formulation**. It evolves both 
%   displacements and element-wise rotation matrices over a sequence of time steps, applying
%   Newton–Raphson iterations at each step to enforce equilibrium.
%
% USAGE:
%   [DATA, CONVERGED, QrotTIME] = NONLINEARstaticLARGE_COROT_LRss(DATA, DISP_CONDITIONS, VAR, ...
%                                                                 OPERFE, SNAP, Fbody, Ftrac, MATPRO)
%
% INPUTS:
%   - DATA            : Structure with general simulation settings, time integration, and solver parameters.
%   - DISP_CONDITIONS : Structure with displacement BCs, projection operators (A, G), and affine deformation vectors.
%   - VAR             : Field variable structure containing displacements, residuals, and internal data.
%   - OPERFE          : Structure containing FE operators, rotation mappings, and element configurations.
%   - SNAP            : Structure for storing snapshot data (displacements, forces, etc.).
%   - Fbody           : Structure with time-dependent body force loadings (U and a).
%   - Ftrac           : Structure with time-dependent traction forces (U and a).
%   - MATPRO          : Structure with material properties (e.g., elasticity, plasticity).
%
% OUTPUTS:
%   - DATA        : Updated DATA structure (time step progression, snapshot info, etc.).
%   - CONVERGED   : Flag indicating if the last time step successfully converged (1 = yes, 0 = no).
%   - QrotTIME    : Cell array storing rotation matrices (`Qrot`) for each time step (used for postprocessing).
%
% FUNCTIONALITY:
%   - Initializes and updates affine-displacement boundary conditions.
%   - Applies external forces, possibly decomposed via SVD (space-time separation).
%   - At each time step, computes the displacement increment and calls:
%       `NewtonRapshonStaticLarge_COROT_LRss` to update both displacements and rotation matrices.
%   - Stores snapshots of the field variables and rotations for later processing or model reduction.
%
% SPECIAL FEATURES:
%   - Designed for corotational formulations in the EIFEM framework.
%   - Maintains per-time-step history of rotation matrices (`Qrot`) for use in visualization or projection.
%   - Compatible with periodic boundary conditions via affine mappings (`A`, `G`).
%
% REFERENCES:
%   - Test problem and example implementation:
%     /TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
%   - Theoretical formulation and explanation:
%     /PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf
%
% AUTHOR:
%   Joaquín A. Hernández, UPC – CIMNE
%   Date: 7-Feb-2025, Honest Greens Pedralbes, Barcelona
%   Comments created by ChatGPT on 12-May-2025
%
% DEPENDENCIES:
%   - NewtonRapshonStaticLarge_COROT_LRss
%   - StoreInfoSnapshots
%   - DefaultField (utility function for safe optional field access)
%
% ---------------------------------------------------------------------------------------------------



% Adaptation of NONLINEARstaticLARGE.m/NONLINEARstaticLARGE_COROT2  
% Small strains/Large rotations (or Small rotations/Large strains)
% JAHO, 7-feb-2025, FRIDAY, 10:48, Honest Greens Pedralbes,  Barcelona.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
% and
%/home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf
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
DATA = DefaultField(DATA,'PROCEED_WITH_NEGATIVE_JACOBIANS',1) ;  

% ROTATION MATRICES, Initialization
% ----------------------------------------
% \DiagC{\QrotALL} = \DiagC{\QrotINIall}
Qrot = OPERFE.QrotINI;  % THIS IS THE MATRIX CONTAINING THE ROTATION MATRICES
% OF ALL THE EIF ELEMENTS. This matrix is constructed in 
%  B_N_matricesEIFEbubCOROT_LRss.m

QrotTIME = cell(1,nsteps) ;   % Here we store all the rotation matrices (for post-process purposes)
QrotTIME{1} = Qrot ;  
% profile on 
while istep <=nsteps
    fprintf(1,'--------------------------------------------------\n')
    fprintf(1,['istep = %d      t_i = %4.5f \n'],[istep,DATA.STEPS(istep)])
    fprintf(1,'--------------------------------------------------\n')
    % 1.a) Initialization displacements    
    % Prescribed macroscopic deformation at  time t_n+1, global coordinates
    dR = DISP_CONDITIONS.dR.U*DISP_CONDITIONS.dR.a(:,istep);
    Delta_dR = dR- DISP_CONDITIONS.dR.U*DISP_CONDITIONS.dR.a(:,istep-1);
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
 
   [VAR,CONVERGED,~,Qrot ]= NewtonRapshonStaticLarge_COROT_LRss(DATA,OPERFE,VAR,MATPRO,DOFl,Qrot,Delta_dR,DOFr) ;
  
    
    QrotTIME{istep} = Qrot ; 
    % Store snaphots
    
    [SNAP,DATA,icluster] = StoreInfoSnapshots(istep,icluster,SNAP,VAR,DATA,CONVERGED) ;
    if CONVERGED == 0  
         
        disp('Convergence error  ....')
        
        
        
        
        pause
        break
    end
    
    istep = istep + 1;
end


%profile report