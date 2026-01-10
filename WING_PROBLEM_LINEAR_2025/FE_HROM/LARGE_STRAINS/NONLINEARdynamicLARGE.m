function [DATA,CONVERGED]=NONLINEARdynamicLARGE(DATA,DISP_CONDITIONS,VAR,OPERFE,SNAP,Fbody,Ftrac,MATPRO)
%%%%
% -------------------------------------------------------------------------
% COMMENTS (generated automatically by ChatGPT on 7-Nov-2025)
%
% PURPOSE:
%   Time-marching solver for large-deformation, possibly nonlinear dynamics.
%   Advances the solution from step 2 to nsteps using a (quasi-)static
%   Newton–Raphson at each time level (with inertia/damping handled via
%   the chosen time-integration operators in VAR). Main tasks per step:
%     1) Build/refresh inertial & damping terms (if ISDYNAMIC==1).
%     2) Enforce prescribed (possibly affine/periodic) displacements.
%     3) Assemble non-follower external forces (body + tractions).
%     4) Solve the nonlinear equilibrium by Newton–Raphson:
%           NewtonRapshonStaticLarge(DATA, OPERFE, VAR, MATPRO, DOFl)
%     5) Store snapshots/iterates and update step counters.
%
% INPUTS:
%   DATA
%     • ISDYNAMIC                 : 1 → dynamic run with inertia/damping operators;
%                                   0 → quasi-static stepping.
%     • InitialAccelerationFromResidual : if 1, initializes acceleration from residual
%                                         via IniAccelFromResidualGET at t = t₁.
%     • CECM_ONLY_FOR_NONLINEAR_STRESSES : hyperreduction toggle for stress terms.
%     • NEWTON_RAPHSON.ENFORCE_BOTH_TOLERANCES_RESIDUAL_AND_DISPLACEMENTS : NR policy flag.
%     • MESH.ndim / MESH.ndimFINE : working dimension (ndim ← ndimFINE).
%     • STEPS                     : vector of time instants {t₁,…,tₙ}; loop starts at istep=2.
%     • TIME_INTEGRATION_PARAMETERS: scheme data used by InertDampExternalForces.
%     • SOLVER_IS_ITERATIVE       : (defaulted to 0 here) passed to the NR solver.
%
%   DISP_CONDITIONS
%     • DOFr, DOFl      : indices of prescribed/free dofs.
%     • DOFm (optional) : master dofs for affine/periodic constraints.
%     • dR.U, dR.a(:,k) : prescribed displacement increments at step k.
%     • A, G (optional) : operators for affine/periodic BCs; applied as
%                         d_R := dR.U*dR.a(:,k) + G*VAR.DISP(DOFm).
%
%   VAR
%     • State variables at the current (converged) time step k−1, including
#       DISP, velocities/accelerations (tilde fields for predictor),
%       and force-like arrays (FINT, FEXT, etc.). Updated in-place.
%
%   OPERFE
%     • Assembly/constraint operators; this routine mirrors DOF maps (DOFr, DOFm)
%       into OPERFE and passes A, G, uBAR (when applicable) to the NR solver.
%
%   SNAP
%     • Snapshot container updated via StoreInfoSnapshots per step/cluster.
%
%   Fbody, Ftrac
%     • Low-rank force representations: F = U*a(:,k). Tractions optional.
%
%   MATPRO
%     • Material/constitutive data consumed by the NR solver.
%
% OUTPUTS:
%   DATA
%     • Updated with per-step info (DeltaT, istep) and any solver metadata
%       added by NewtonRapshonStaticLarge and snapshot utilities.
%   CONVERGED
%     • 1 if the last processed step converged; 0 if a convergence failure
%#       occurred (loop is broken immediately on failure).
%
% STEP-BY-STEP FLOW:
%   0) Optional initial acceleration from residual (t = t₁) if flags demand.
%   1) For istep = 2..nsteps:
%        - DeltaT ← t_k − t_{k−1}
%        - UpdateIterativeSnapshot(DATA, VAR, iterk=1) for logging the previous converged state
%        - If ISDYNAMIC:
%            VAR ← InertDampExternalForces(VAR, TIME_INTEGRATION_PARAMETERS, DeltaT, OPERFE)
%            VAR.DISP(DOFl) ← VAR.DISPtilde(DOFl)     % predictor on free dofs
%        - Prescribed displacements at t_k:
%            dR := dR.U*dR.a(:,k); if DOFm present → dR += G*VAR.DISP(DOFm)
%            VAR.DISP(DOFr) ← dR;     OPERFE.uBAR ← dR
%        - External forces (non-follower):
%            VAR.FEXT ← Fbody.U*a_body(:,k) [+ Ftrac.U*a_trac(:,k) if available]
%        - Nonlinear solve:
%            DATA.DeltaT ← DeltaT; DATA.istep ← istep
%            [VAR, CONVERGED, DATA] ← NewtonRapshonStaticLarge(DATA, OPERFE, VAR, MATPRO, DOFl)
%        - StoreInfoSnapshots(istep, icluster, SNAP, VAR, DATA, CONVERGED)
%        - If ~CONVERGED → print diagnostic and break
%        - istep ← istep + 1
%
% ASSUMPTIONS / SANITY:
%   - STEPS is strictly increasing; istep starts at 2 (state at step 1 assumed known).
%   - Dimensions of F*.U and a(:,k) conform with VAR.DISP length (ndof).
%   - DOFr, DOFl, DOFm are consistent and non-overlapping.
%   - Affine/periodic constraints are enforced through A, G, and uBAR in OPERFE;
%     the exact constraint algebra is handled by the NR solver.
%
% NUMERICAL NOTES:
%   - Convergence criteria and line-search/damping are encapsulated in
%     NewtonRapshonStaticLarge (optionally enforcing residual & displacement norms).
%   - Hyperreduction flags (e.g., CECM_ONLY_FOR_NONLINEAR_STRESSES) are passed
%     via DATA and should be honored inside the assembly routines called by NR.
%   - Snapshot collection is decoupled (StoreInfoSnapshots) to keep the loop lean.
%
% REFERENCES (paths in developer workspace):
%   - See comments in code for pointers to test notebooks (e.g., 19_ExactLinearStiff.mlx)
%     and periodic/affine BC examples.
%
% AUTHOR / HISTORY:
%   Comments clarification: 7-Nov-2025
%   JAHO — Joaquín A. Hernández — jhortega@cimne.upc.edu
% -------------------------------------------------------------------------

if nargin == 0
    load('tmp3.mat')
end

% Initial acceleration
% -------------------
DATA = DefaultField(DATA,'InitialAccelerationFromResidual',0) ;
DATA = DefaultField(DATA,'CECM_ONLY_FOR_NONLINEAR_STRESSES',0) ; % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx

DATA.NEWTON_RAPHSON = DefaultField(DATA.NEWTON_RAPHSON,'ENFORCE_BOTH_TOLERANCES_RESIDUAL_AND_DISPLACEMENTS',0) ; % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx

DATA.MESH = DefaultField(DATA.MESH,'ndimFINE',DATA.MESH.ndim) ;
DATA.MESH.ndim = DATA.MESH.ndimFINE;

DOFr = DISP_CONDITIONS.DOFr ;
DOFm = [] ; 
if isfield(DISP_CONDITIONS,'DOFm')
        DOFm = DISP_CONDITIONS.DOFm ;
    OPERFE.A = DISP_CONDITIONS.A;
    OPERFE.G = DISP_CONDITIONS.G;
    OPERFE.DOFr =   DOFr ;
    OPERFE.DOFm =DOFm;
else
    OPERFE.DOFm = [] ;
    OPERFE.A = [];
end


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
    end
    
    
    
    
    
    % 1.b) Nodal  External forces (non-follower)
%    VAR.FEXT = Fbody.U*Fbody.a(:,istep) + Ftrac.U*Ftrac.a(:,istep) ;
%if ~isempty(Fbody.U)
    VAR.FEXT =   Fbody.U*Fbody.a(:,istep) ;
%end
if ~isempty(Ftrac.U)
    VAR.FEXT = VAR.FEXT +  Ftrac.U*Ftrac.a(:,istep) ;
end
    
    % *********************************
    %%%% NEWTON-RAPHSON ALGORITHM
    % *****************************
    
    DATA.DeltaT = DeltaT ; DATA.istep = istep  ;
    [VAR,CONVERGED,DATA ]= NewtonRapshonStaticLarge(DATA,OPERFE,VAR,MATPRO,DOFl) ;
    
    % Store snaphots
    [SNAP,DATA,icluster] = StoreInfoSnapshots(istep,icluster,SNAP,VAR,DATA,CONVERGED) ;
    if CONVERGED == 0
        disp('Convergence error ....')
        break
    end
    istep = istep + 1;
end