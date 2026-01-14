function K = GetStiffnessMatrix_COROT(OPERFE,DATA,VAR,FgradST,celastST,KcLINq,BstQ,D_QrotALL)
% Adaptation of GetStiffnessMatrix.m to the co-rotational method
% JAHO, 29-Oct-2024, Tuesday, Honest Greens, Pedralbes Center, Barcelona
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
% ----------------------------------------------------------------------

% if  ~isempty(OPERFE.KinternalFORCES_given)
%     % Stiffness matrix given by the user
%     % This was implemented when testing the Dynamic Mode Decomposition 
%     % Decomposition, see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/DynamicModeDecomposition/01_ForcedVibration_2D/RunDMD_aux.mlx
%     if ~isstruct(OPERFE.KinternalFORCES_given)
%         K = OPERFE.KinternalFORCES_given ;
%     else
%         iMATRIX = OPERFE.KinternalFORCES_given.INTERVALS(DATA.istep) ;
%         K = OPERFE.KinternalFORCES_given.MATRICES{iMATRIX} ;
%     end
% else
    if DATA.SMALL_STRAIN_KINEMATICS == 0
        K = KstiffLargeStrains_COROT(OPERFE,VAR.PK2STRESS,FgradST,DATA.MESH.ndim,celastST,DATA,KcLINq,BstQ) ;
    else
        error('Option not implemented')
        K = KstiffSmallStrains(OPERFE,FgradST,DATA.MESH.ndim,celastST,DATA) ;
    end
% end