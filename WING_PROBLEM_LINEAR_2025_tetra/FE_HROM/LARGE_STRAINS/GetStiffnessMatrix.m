function K = GetStiffnessMatrix(OPERFE,DATA,VAR,FgradST,celastST)
%--------------------------------------------------------------------------
% K = GetStiffnessMatrix(OPERFE, DATA, VAR, FgradST, celastST)
%
% PURPOSE:
%   Computes or retrieves the global stiffness matrix K based on the current
%   kinematic configuration, material tangent, and user options. This is a
%   wrapper function that selects between:
%     - a user-provided stiffness matrix (for testing, e.g., DMD),
%     - the large strain stiffness matrix,
%     - the small strain stiffness matrix.
%
% DESCRIPTION:
%   This function is used during the solution procedure to compute internal
%   forces and update the tangent matrix within the Newton-Raphson scheme.
%   It automatically selects the appropriate stiffness formulation depending
%   on the problem type (small vs large strain) and user-provided options.
%
% INPUT:
%   OPERFE     : Structure containing finite element operators and, optionally:
%                - KinternalFORCES_given : either a constant K matrix or a struct
%                  with time-dependent matrices defined by INTERVALS and MATRICES
%   DATA       : Structure with analysis settings, including:
%                - SMALL_STRAIN_KINEMATICS : boolean flag
%                - istep : current time/load step (for INTERVALS usage)
%   VAR        : Structure with current solution variables
%                - PK2STRESS : second Piola–Kirchhoff stresses (used in large strains)
%   FgradST    : Deformation gradient at integration points
%   celastST   : Consistent tangent operator (Voigt format) at Gauss points
%
% OUTPUT:
%   K          : Assembled stiffness matrix for the current configuration
%
% NOTES:
%   - If OPERFE.KinternalFORCES_given is set, it overrides the automatic assembly.
%   - In large strain problems, PK2 stresses are used in geometric stiffness.
%   - In small strain problems, linearized elasticity is assumed.
%
% RELATED FILES:
%   - KstiffLargeStrains.m : Large strain stiffness matrix computation
%   - KstiffSmallStrains.m : Small strain stiffness matrix computation
%   - RunDMD_aux.mlx       : Context for DMD testing with prescribed stiffness
%--------------------------------------------------------------------------

if  ~isempty(OPERFE.KinternalFORCES_given)
    % Stiffness matrix given by the user
    % This was implemented when testing the Dynamic Mode Decomposition 
    % Decomposition, see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/DynamicModeDecomposition/01_ForcedVibration_2D/RunDMD_aux.mlx
    if ~isstruct(OPERFE.KinternalFORCES_given)
        K = OPERFE.KinternalFORCES_given ;
    else
        iMATRIX = OPERFE.KinternalFORCES_given.INTERVALS(DATA.istep) ;
        K = OPERFE.KinternalFORCES_given.MATRICES{iMATRIX} ;
    end
else
    if DATA.SMALL_STRAIN_KINEMATICS == 0
        K = KstiffLargeStrains(OPERFE,VAR.PK2STRESS,FgradST,DATA.MESH.ndim,celastST,DATA) ;
    else
        K = KstiffSmallStrains(OPERFE,FgradST,DATA.MESH.ndim,celastST,DATA) ;
    end
end