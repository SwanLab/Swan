function INICOND = InitialDisplacementsDueToGravity(DATA,INICOND,DISP_CONDITIONS,Fbody,Ftrac,OPERFE,MATPRO)
%%%%
% -------------------------------------------------------------------------
% COMMENTS (generated automatically by ChatGPT on 7-Nov-2025)
%
% PURPOSE:
%   Compute the initial displacement field induced by gravity (and other
%   time-step-1 loads) by solving a static equilibrium problem via a
%   Newton–Raphson procedure. The routine:
%     1) Initializes state variables (INITIALIZATIONvar).
%     2) Applies prescribed displacement increments at step 1.
%     3) Builds external forces from body loads and tractions at step 1.
%     4) Solves a static nonlinear problem (large-displacement kinematics)
%        to obtain VAR.DISP, which is then stored in INICOND.DISP.
%
% INPUTS:
%   DATA             - global settings; relevant flags include:
%                        • DATA.INTERNAL_FORCES_USING_precomputed_Kstiff (bool)
%                        • (internally set) DATA.SKIP_PART_STORE = 1
%   INICOND          - structure carrying initial conditions (updated here).
%   DISP_CONDITIONS  - structure with boundary conditions:
%                        • dR.U, dR.a(:,1) → prescribed dof increments at step 1
%                        • DOFr, DOFl      → free/fixed dof index sets
%   Fbody, Ftrac     - low-rank representations of loads:
%                        • F*.U, F*.a(:,1) so that F = U * a(:,1)
%   OPERFE           - operators/flags for FE assembly:
%                        • KinternalFORCES_given (optional precomputed data)
%                        • DOFm (optional)
%   MATPRO           - material/constitutive data used by the solver.
%
% OUTPUTS:
%   INICOND          - updated with:
%                        • INICOND.DISP = VAR.DISP (equilibrium under gravity at step 1)
%
% METHOD / STEPS:
%   - Initialization:
%       DATA = DefaultField(DATA,'INTERNAL_FORCES_USING_precomputed_Kstiff',0);
%       [DATA,VAR,~] = INITIALIZATIONvar(DATA,INICOND);
%   - Apply prescribed displacements (if any) at step 1:
%       VAR.DISP(DISP_CONDITIONS.DOFr) ← dR.U * dR.a(:,1)
%   - External forces at step 1:
%       VAR.FEXT = Fbody.U*Fbody.a(:,1) + Ftrac.U*Ftrac.a(:,1)
%   - Static solve (large-displacement, no time integration):
%       DATAstatic ← DATA; set ISDYNAMIC=0, SOLVER_IS_ITERATIVE=0
%       [VAR, CONVERGED] = NewtonRapshonStaticLarge(DATAstatic, OPERFE, VAR, MATPRO, DISP_CONDITIONS.DOFl)
%   - Set initial condition:
%       INICOND.DISP ← VAR.DISP
%
% ASSUMPTIONS / SANITY:
%   - Step index fixed to istep = 1 (initial gravitational load state).
%   - Dimensions of Fbody.U/Ftrac.U and coefficients a(:,1) are conformable
%     with VAR.DISP (ndof × 1).
%   - DOFr/DOFl are consistent, non-overlapping index sets.
%   - Convergence handling is delegated to NewtonRapshonStaticLarge; no
%     explicit fallback is implemented here.
%
% AUTHOR / HISTORY:
%   Comments clarification: 7-Nov-2025
%   JAHO — Joaquín A. Hernández — jhortega@cimne.upc.edu
% -------------------------------------------------------------------------

if nargin == 0
    load('tmp.mat')
end

DATA.SKIP_PART_STORE = 1;
DATA = DefaultField(DATA,'INTERNAL_FORCES_USING_precomputed_Kstiff',0)  ;

[DATA,VAR,~] = INITIALIZATIONvar(DATA,INICOND)   ;

istep = 1;

if ~isempty(DISP_CONDITIONS.dR.U)
    dR = DISP_CONDITIONS.dR.U*DISP_CONDITIONS.dR.a(:,istep);
    VAR.DISP(DISP_CONDITIONS.DOFr) = dR;
end

% 1.b) External forces
VAR.FEXT = Fbody.U*Fbody.a(:,istep) + Ftrac.U*Ftrac.a(:,istep) ;

OPERFE = DefaultField(OPERFE,'KinternalFORCES_given',[]) ;
OPERFE = DefaultField(OPERFE,'DOFm',[]) ;


% *********************************
%%% NEWTON-RAPHSON ALGORITHM
% *****************************+
DATAstatic = DATA ;
DATAstatic.ISDYNAMIC = 0 ;
DATAstatic.SOLVER_IS_ITERATIVE = 0 ;
DATAstatic.SNAP_ITER = [] ;
[VAR,CONVERGED ]= NewtonRapshonStaticLarge(DATAstatic,OPERFE,VAR,MATPRO,DISP_CONDITIONS.DOFl) ;

% Re-setting initial condition for displacements
INICOND.DISP = VAR.DISP ;