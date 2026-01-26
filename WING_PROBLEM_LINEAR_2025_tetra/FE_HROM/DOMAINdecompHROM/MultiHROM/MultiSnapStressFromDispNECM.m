function [SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,VAR,ReactFORCE,SNAPstressSTWOproj_LOC_INELAST,SNAPstressSTWOproj_LOC_ELAST] = ...
    MultiSnapStressFromDispNECM(VAR,d,DATA,VARint_n,OTHER_output,OPERFE,MATPRO,iloc,...
    SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,iproj,INFO_RVE,Fbody_loc,Ftrac_loc,SNAPstressSTWOproj_LOC_INELAST,...
    SNAPstressSTWOproj_LOC_ELAST)
%--------------------------------------------------------------------------
% [SNAPstressPonePROJ_LOC, SNAPstressSTWOproj_LOC, VAR, ReactFORCE, ...
%  SNAPstressSTWOproj_LOC_INELAST, SNAPstressSTWOproj_LOC_ELAST] = ...
%   MultiSnapStressFromDispNECM(VAR, d, DATA, VARint_n, OTHER_output, OPERFE, ...
%                               MATPRO, iloc, SNAPstressPonePROJ_LOC, ...
%                               SNAPstressSTWOproj_LOC, iproj, INFO_RVE, ...
%                               Fbody_loc, Ftrac_loc, ...
%                               SNAPstressSTWOproj_LOC_INELAST, ...
%                               SNAPstressSTWOproj_LOC_ELAST)
%
% PURPOSE:
%   Computes stress snapshots (PK1 and PK2), elastic and inelastic stress
%   contributions, and reaction forces at a given location (iloc) from a
%   prescribed displacement vector. This routine extends
%   SnapStressFromDispLOC by separating elastic and inelastic stresses,
%   which is required for nonlinear constitutive laws in reduced-order
%   modeling (e.g., EIFEM, HROM).
%
% FUNCTIONALITY:
%   - Updates the displacement field VAR.DISP with input d.
%   - Initializes or updates internal variables (plastic strains, yield
%     stress, etc.) depending on location index:
%       * iloc = 1 → load initial conditions from OTHER_output.INICOND.
%         Specific indexing into Gauss points is performed via INFO_RVE.
%       * iloc > 1 → reuse stored internal variables in VAR.
%   - Computes elastic stresses:
%       * Builds VAR_for_disp_zero with zero displacements (virgin state).
%       * Calls StressesFromDisplacementsVAR to extract the elastic tangent
%         stiffness celastST_initial.
%       * Forms elastic stresses as STRESS_ELASTIC = C_elast * strain(d).
%   - Computes stresses with material history:
%       * Calls StressesFromDisplacementsVARincre to update VAR with full
%         constitutive response (PK2, PK1).
%   - Updates stress snapshots:
%       * SNAPstressSTWOproj_LOC{iloc}       → total PK2 stresses.
%       * SNAPstressSTWOproj_LOC_ELAST{iloc} → purely elastic stresses.
%       * SNAPstressSTWOproj_LOC_INELAST{iloc} → inelastic stresses 
%           (PK2 - elastic part), unless TRAIN_WITH_TOTAL_STRESSES=1, 
%           in which case it stores total PK2.
%       * SNAPstressPonePROJ_LOC{iproj}      → PK1 stresses if available,
%           otherwise PK2.
%   - Internal forces:
%       * VAR.FINT = InternalForces(...) from PK stresses.
%   - External forces:
%       * Computed from body forces (Fbody_loc) and tractions (Ftrac_loc).
%   - Reaction forces:
%       * ReactFORCE = VAR.FINT – Fext.
%
% INPUTS:
%   VAR        : State variables structure (displacements, stresses, internals).
%   d          : Displacement vector at current configuration.
%   DATA       : Global analysis settings.
%   VARint_n   : Internal variables at previous step/iloc.
%   OTHER_output : Structure with initial condition fields.
%   OPERFE     : Finite element operators (strain-displacement matrices, etc.).
%   MATPRO     : Material property structure (constitutive law data).
%   iloc       : Local snapshot index (location).
%   SNAPstressPonePROJ_LOC : Cell array storing PK1 stress snapshots.
%   SNAPstressSTWOproj_LOC : Cell array storing PK2 stress snapshots.
%   iproj      : Projection index (can differ from iloc).
%   INFO_RVE   : Microstructural/RVE information (Gauss indices, DOFs).
%   Fbody_loc  : Body force contribution structure.
%   Ftrac_loc  : Traction contribution structure.
%   SNAPstressSTWOproj_LOC_INELAST : Cell array to store inelastic stress snapshots.
%   SNAPstressSTWOproj_LOC_ELAST   : Cell array to store elastic stress snapshots.
%
% OUTPUTS:
%   SNAPstressPonePROJ_LOC         : Updated with PK1 (or PK2 if PK1 empty).
%   SNAPstressSTWOproj_LOC         : Updated with total PK2 stresses.
%   VAR                            : Updated state variables with stresses.
%   ReactFORCE                     : Reaction forces = internal – external.
%   SNAPstressSTWOproj_LOC_INELAST : Inelastic stresses (PK2 – elastic).
%   SNAPstressSTWOproj_LOC_ELAST   : Elastic stresses.
%
% CONTEXT:
%   This function is used in nonlinear ECM/HROM/EIFEM workflows to build
%   stress snapshots, separating elastic and inelastic parts, for
%   hyperreduction and projection-based model order reduction. Reaction
%   forces are computed for consistency checks in static problems.
%--------------------------------------------------------------------------
% JAHO, comments by ChatGPT-5, 23th August 2025

% Copy of MultiSnapStressFromDispLOC.m
% Adapted for extracting nonlinear stresses
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx

if nargin == 0
    load('tmp.mat')
end



VAR.DISP  =d ;


if isempty(DATA.ListFieldInternalVariables)
    VARint_n = [] ;
else
    if iloc == 1
        for ivar = 1:length(DATA.ListFieldInternalVariables)
            FLOC = DATA.ListFieldInternalVariables{ivar};
            
            
            switch FLOC
                % Temporal amendment !!
                case {'YieldStress','InternalVarStrain'}
                    VARint_n.(FLOC) = OTHER_output.INICOND.(FLOC)(INFO_RVE.GaussINDEX_scalarSTR,:)  ;
                    
                case {'PlasticStrains'}
                    VARint_n.(FLOC) = OTHER_output.INICOND.(FLOC)(INFO_RVE.GaussINDEX_stress,:)  ;
            end
            VAR.(FLOC) =  VARint_n.(FLOC)  ;
            
            
        end
    else
        for ivar = 1:length(DATA.ListFieldInternalVariables)
            FLOC = DATA.ListFieldInternalVariables{ivar};
            VARint_n.(FLOC) = VAR.(FLOC)  ;
        end
    end
    
    
end


VAR_for_disp_zero = VAR;
VAR_for_disp_zero.DISP =zeros(DATA.MESH.ndof,1) ;  % Global elasticity matrix for disp = zero (virgin material)
DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES = 0;
[~,celastST_initial,~,~] = StressesFromDisplacementsVAR(OPERFE,VAR_for_disp_zero,MATPRO,DATA,VARint_n) ;

STRAIN = OPERFE.Bst*d;
celastST_initial = ConvertBlockDiag(celastST_initial) ;
STRESS_ELASTIC = (celastST_initial*STRAIN);


% STRESSES FROM HISTORY OF DISPLACEMENTS
[VAR,~,~,~] = StressesFromDisplacementsVARincre(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
SNAPstressSTWOproj_LOC{iloc} = VAR.PK2STRESS ;

TRAIN_WITH_TOTAL_STRESSES = 0; % sEE /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx

SNAPstressSTWOproj_LOC_ELAST{iloc} = STRESS_ELASTIC ;

if TRAIN_WITH_TOTAL_STRESSES == 1
    SNAPstressSTWOproj_LOC_INELAST{iloc} = VAR.PK2STRESS  ;  % This was implemented when
    % developing the approach (to see what it failed in first place), see
    % sEE /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
else
    SNAPstressSTWOproj_LOC_INELAST{iloc} = VAR.PK2STRESS   -STRESS_ELASTIC ;  % This is the efficient option
end
if  isempty(VAR.PK1STRESS)
    SNAPstressPonePROJ_LOC{iproj} = VAR.PK2STRESS ; ;
else
    SNAPstressPonePROJ_LOC{iproj} = VAR.PK1STRESS ;
end


% Internal forces (FOR ALL TIME STEPS)
% ---------------
VAR.FINT = InternalForces(OPERFE,VAR.PK1STRESS,VAR.PK2STRESS,DATA) ;

% ------------------------------------------------------------------------------
% REACTIVE FORCES (DIFFERENCE BETWEEN INTERNAL FORCES AND EXTERNAL FORCES)
% ------------------------------------------------------------------------------
% ONLY STATIC CASE SO FAR (8-Feb-2023)
% --------------------------------------------------
Fext = 0 ;
if ~isempty(Fbody_loc.U)
    Fext = Fext + Fbody_loc.U(INFO_RVE.DOFS_globNUM,:)*Fbody_loc.a;
end
if ~isempty(Ftrac_loc.U)
    Fext = Fext + Ftrac_loc.U(INFO_RVE.DOFS_globNUM,:)*Ftrac_loc.a;
end

ReactFORCE= VAR.FINT-Fext ;

