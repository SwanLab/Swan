function [SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,VAR,ReactFORCE,SNAPstressPONEproj_LOC_NONL,SNAPstressPONEproj_LOC_LINEAR] = ...
    MultiSnapStressFromDispNECMlarg(VAR,d,DATA,VARint_n,OTHER_output,OPERFE,MATPRO,iloc,...
    SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,iproj,INFO_RVE,Fbody_loc,Ftrac_loc,SNAPstressPONEproj_LOC_NONL,...
    SNAPstressPONEproj_LOC_LINEAR)
% ---------------------------------------------------------------------------------------------------
% FUNCTION: MultiSnapStressFromDispNECMlarg
% ---------------------------------------------------------------------------------------------------
% PURPOSE:
%   Given a reduced (projected) displacement field, this function computes:
%     - First Piola–Kirchhoff (PK1) stresses (linear and nonlinear contributions)
%     - Second Piola–Kirchhoff (PK2) stresses
%     - Internal forces and reactive forces (Fint - Fext)
%
%   It is designed for use in **Nonlinear Empirical Cubature Methods (NECM)** and **CECM-style**
%   stress projection pipelines under **large strain kinematics**, with the goal of distinguishing
%   between linear (elastic) and nonlinear (inelastic) components of stress fields.
%
% USAGE:
%   [SNAPstressPonePROJ_LOC, SNAPstressSTWOproj_LOC, VAR, ReactFORCE, ...
%    SNAPstressPONEproj_LOC_NONL, SNAPstressPONEproj_LOC_LINEAR] = ...
%       MultiSnapStressFromDispNECMlarg(VAR, d, DATA, VARint_n, OTHER_output, OPERFE, MATPRO, ...
%                                       iloc, SNAPstressPonePROJ_LOC, SNAPstressSTWOproj_LOC, ...
%                                       iproj, INFO_RVE, Fbody_loc, Ftrac_loc, ...
%                                       SNAPstressPONEproj_LOC_NONL, SNAPstressPONEproj_LOC_LINEAR)
%
% INPUTS:
%   - VAR                         : Structure containing all current field variables (e.g., displacements).
%   - d                           : Projected displacement vector.
%   - DATA                        : Main simulation parameters, flags, and configuration.
%   - VARint_n                    : Structure with internal variable history (used in elasto-plasticity).
%   - OTHER_output                : Structure with problem-specific precomputed variables (e.g., INICOND).
%   - OPERFE                      : Finite element operators and mappings (e.g., shape functions, Booleans).
%   - MATPRO                      : Material properties (e.g., elastic moduli, constitutive models).
%   - iloc                        : Local index for the current snapshot within a case.
%   - SNAPstressPonePROJ_LOC      : Cell to store PK1 stress snapshots (projected).
%   - SNAPstressSTWOproj_LOC      : Cell to store PK2 stress snapshots (projected).
%   - iproj                       : Index of the current projection case.
%   - INFO_RVE                    : Structure with indices for selected DOFs and Gauss points.
%   - Fbody_loc, Ftrac_loc        : Localized body and traction force vectors for computing Fext.
%   - SNAPstressPONEproj_LOC_NONL : Cell to store inelastic (nonlinear) PK1 stress component.
%   - SNAPstressPONEproj_LOC_LINEAR : Cell to store linear elastic PK1 stress component.
%
% OUTPUTS:
%   - SNAPstressPonePROJ_LOC      : Updated PK1 stress snapshots (total).
%   - SNAPstressSTWOproj_LOC      : Updated PK2 stress snapshots (total).
%   - VAR                         : Updated field variables including stress fields.
%   - ReactFORCE                  : Residual forces (Fint - Fext) at selected DOFs.
%   - SNAPstressPONEproj_LOC_NONL : Snapshot of nonlinear component of PK1 stresses.
%   - SNAPstressPONEproj_LOC_LINEAR : Snapshot of linear (elastic) component of PK1 stresses.
%
% FUNCTIONALITY:
%   - Forces small strain assumption to compute a consistent linear elastic reference stress field.
%   - Computes full nonlinear PK1 and PK2 stresses using large strain kinematics.
%   - Computes internal forces from stresses.
%   - Computes external forces from provided body/traction data.
%   - Outputs stress differences for separation of elastic/inelastic response (for training or error analysis).
%
% NOTES:
%   - Plasticity with history-dependent variables is **not yet supported** under large strain (see TODO).
%   - Internal variables are used only in small strain approximation for elastic reference.
%   - Compatible with `GetStressesAndReactForces_bubNECMlarg` for stress basis projection.
%
% REFERENCES:
%   - Development and theory in:
%     /TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
%   - Related conceptual discussion in:
%     /PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf
%
% AUTHOR:
%   Joaquín A. Hernández, UPC – CIMNE  
%   Date: 23-May-2024, Campus Nord, Barcelona
%   Comments by ChatGPT-4, 12-May-2025
%
% DEPENDENCIES:
%   - StressesFromDisplacementsVAR
%   - StressesFromDisplacementsVARincre
%   - InternalForces
%   - DefaultField
%
% ---------------------------------------------------------------------------------------------------




% Copy of MultiSnapStressFromDispNECM.m
% Adapted for extracting nonlinear stresses in large strains 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx

if nargin == 0
    load('tmp3.mat')
end

VAR.DISP  =d ;

if isempty(DATA.ListFieldInternalVariables)
    VARint_n = [] ;
else
    error('Option not ready for plasticity + large strains')
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

SMALL_STRAIN_KINEMATICS = DATA.SMALL_STRAIN_KINEMATICS  ; 
%VAR_for_disp_zero = VAR;
%VAR_for_disp_zero.DISP =zeros(DATA.MESH.ndof,1) ;  % Global elasticity matrix for disp = zero (virgin material)
DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES = 0;
DATA.SMALL_STRAIN_KINEMATICS = 1; % We force small strains kinematics (linear)
DATA.NO_USE_Deformation_gradient_in_Small_Strains = 0 ; 
[VAR_OUT,~,~,~] = StressesFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
 
STRESS_ELASTIC = VAR_OUT.PK1STRESS ; % These are the linear PK1 stresses (4 or 9 entries per integration point)


% STRESSES FROM HISTORY OF DISPLACEMENTS
% Our goal is to compare PK2 stresses, but ECM is calculated with PK1, so
% we need both 
DATA.SMALL_STRAIN_KINEMATICS = SMALL_STRAIN_KINEMATICS;
[VAR,~,~,~] = StressesFromDisplacementsVARincre(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
SNAPstressSTWOproj_LOC{iloc} = VAR.PK2STRESS ;

TRAIN_WITH_TOTAL_STRESSES = 0; % sEE /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx

SNAPstressPONEproj_LOC_LINEAR{iloc} = STRESS_ELASTIC ;

if TRAIN_WITH_TOTAL_STRESSES == 1
    SNAPstressPONEproj_LOC_NONL{iloc} = VAR.PK2STRESS  ;  % This was implemented when
    % developing the approach (to see what it failed in first place), see
    % sEE /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
else
   % SNAPstressPONEproj_LOC_NONL{iloc} = VAR.PK2STRESS   -STRESS_ELASTIC ;  % This is the efficient option
    SNAPstressPONEproj_LOC_NONL{iloc} = VAR.PK1STRESS   -STRESS_ELASTIC ;  % This is for large strains 
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

COMPUTE_REACTIONS = 0; 

if COMPUTE_REACTIONS == 1
Fext = 0 ;
if ~isempty(Fbody_loc.U)
    Fext = Fext + Fbody_loc.U(INFO_RVE.DOFS_globNUM,:)*Fbody_loc.a;
end
if ~isempty(Ftrac_loc.U)
    Fext = Fext + Ftrac_loc.U(INFO_RVE.DOFS_globNUM,:)*Ftrac_loc.a;
end

ReactFORCE= VAR.FINT-Fext ;

else
    % Disabled 26-May-2025
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/13_auxALEX.mlx
    ReactFORCE = [] ; 
end