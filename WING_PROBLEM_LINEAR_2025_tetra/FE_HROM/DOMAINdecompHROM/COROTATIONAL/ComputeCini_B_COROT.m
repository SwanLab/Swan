function [ OPERFE]= ComputeCini_B_COROT(DATA,OPERFE,MATPRO,VAR)
%--------------------------------------------------------------------------
% function [OPERFE] = ComputeCini_B_COROT(DATA, OPERFE, MATPRO, VAR)
%
% PURPOSE:
%   Computes and stores the linear elastic constitutive matrix at the
%   initial configuration, required for the Empirical Cubature-based
%   approximation of nonlinear internal forces:
%
%       P = P_linear + P_nonlinear,
%   where only the nonlinear part is approximated using ECM. This setup is 
%   compatible with the corotational approach and the EIFEM method for 
%   problems involving large rotations with small or moderate strains.
%
% USAGE CONTEXT:
%   Used when DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 1. This allows
%   separation of elastic (linear) contributions from the total internal forces.
%
% BEHAVIOR:
%   - In nonlinear (finite strain) mode, the routine switches temporarily 
%     to SMALL_STRAIN_KINEMATICS = 1 to extract the linear tangent stiffness.
%   - The result is stored in OPERFE.D_CtangFlinST as a **diagonal block matrix**.
%   - The diagonalization facilitates efficient evaluation of linear internal forces.
%
% INPUT:
%   DATA      : Global data structure including simulation flags, mesh, and settings.
%   OPERFE    : Finite element operators (contains Bst, BstW, etc.).
%   MATPRO    : Material property structure (stiffness tensors, models).
%   VAR       : Structure with current state variables (displacement, internal vars).
%
% OUTPUT:
%   OPERFE    : Updated with the field `D_CtangFlinST`, which is the global
%               diagonalized elastic tangent matrix for all Gauss points.
%
% REMARKS:
%   - If `DATA.SMALL_STRAIN_KINEMATICS == 1`, the function raises an error,
%     as corotational splitting for purely small strains is not implemented.
%   - This function is an adaptation of `ComputeCini.m` to the context of
%     large rotations (COROT), where a linear reference configuration is needed.
%
% REFERENCES:
%   - See `/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx`
%   - See also related developments in:
%       `/108_EIFEM_metamat/04_EIFEM_NONL.mlx`
%       `/FE_HROM/LARGE_STRAINS/KINEMATICS/KstiffLargeStrains.m`
%
% AUTHOR:
%   Joaquín A. Hernández Ortega, UPC-CIMNE, Balmes 185, Barcelona
%   Created: 28-Oct-2024
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------




% Adaptation of ComputeCini.m to cope with the corotational approach 
% Linear constitutive matrix for all ECM points and elements 
% \DiagC{\CtangFlinST}  \defeq  \diagOL{\CtangFlinSTe{1}}{\CtangFlinSTe{2}}{\cdots}{\CtangFlinSTe{\nelemC}}
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
% JAHO, 28-Oct-2024, Balmes 185, Barcelona 

% ----------------------------------------------
if nargin == 0
    load('tmp1.mat')
end



if DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 1
    
    
    % Internal variables (value at previous time step)
    % ----------------------
    VARint_n = [];
    if ~isempty(DATA.ListFieldInternalVariables)
        for iINTVAR = 1:length(DATA.ListFieldInternalVariables)
            NameIntVar= DATA.ListFieldInternalVariables{iINTVAR} ;
            VARint_n.(NameIntVar) = VAR.(NameIntVar) ;
        end
    end
    SMALL_STRAIN_KINEMATICS = DATA.SMALL_STRAIN_KINEMATICS ;  % =
    
    if SMALL_STRAIN_KINEMATICS == 1
        error('Option not implemented')
         DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES = 0 ;
        [~, OPERFE.celastST_ini,~,~] = StressesFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
        
    else
        % See  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
        % We are interested in the linear part of the deformation
        % This is obtained by setting
        DATA.SMALL_STRAIN_KINEMATICS = 1;
        DATA.NO_USE_Deformation_gradient_in_Small_Strains = 0 ;
        DATA.MESH.ndim = DATA.MESH.ndimFINE ;
        DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES = 0 ;
        [~, celastST_ini,FgradST,~] = StressesFromDisplLIN(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
        
        % In contrast to the small-strains case, here we will turn
        % celastST_ini into a diagonal, block matrix, see more details in
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
        % See also /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/KINEMATICS/KstiffLargeStrains.m
        
        OPERFE.D_CtangFlinST = CelasLARGEmat_allgauss(celastST_ini,FgradST,DATA.MESH.ndim) ;
        
         OPERFE.D_CtangFlinST = ConvertBlockDiag(OPERFE.D_CtangFlinST) ; % Diagonal block matrix
        
        
        
    end
    
    
end



%