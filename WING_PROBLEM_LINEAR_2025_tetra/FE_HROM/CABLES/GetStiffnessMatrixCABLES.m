function K = GetStiffnessMatrixCABLES(OPERFE,DATA,VAR,FgradST,celastST)


 nF = size(celastST,2) ; 

for icomp = 1:nF
    icol = icomp:nF:size(celastST,1) ;
    celastST(icol,:) = bsxfun(@times,celastST(icol,:),OPERFE.wSTs) ; 
end
 
celastST = ConvertBlockDiag(celastST) ; % Diagonal block matrix 

K = OPERFE.Bst'*(celastST*OPERFE.Bst); 



% if  ~isempty(OPERFE.KinternalFORCES_given)
%     % Stiffness matrix given by the user
%     % This was implemented when testing the Dynamic Mode
%     % Decomposition, see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/DynamicModeDecomposition/01_ForcedVibration_2D/RunDMD_aux.mlx
%     if ~isstruct(OPERFE.KinternalFORCES_given)
%         K = OPERFE.KinternalFORCES_given ;
%     else
%         iMATRIX = OPERFE.KinternalFORCES_given.INTERVALS(DATA.istep) ;
%         K = OPERFE.KinternalFORCES_given.MATRICES{iMATRIX} ;
%     end
% else
%     if DATA.SMALL_STRAIN_KINEMATICS == 0
    %   K = KstiffLargeStrains(OPERFE,VAR.PK2STRESS,FgradST,DATA.MESH.ndim,celastST) ;
%     else
%         K = KstiffSmallStrains(OPERFE,FgradST,DATA.MESH.ndim,celastST) ;
%     end
% end