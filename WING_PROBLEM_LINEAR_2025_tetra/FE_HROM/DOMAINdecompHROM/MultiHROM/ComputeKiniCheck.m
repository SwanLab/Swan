function [ OTHER_output,OPERFE]= ComputeKiniCheck(DATA,PARAMETERS,OTHER_output,OPERFE,MATPRO)
%--------------------------------------------------------------------------
% function [OTHER_output, OPERFE] = ComputeKiniCheck(DATA, PARAMETERS, OTHER_output, OPERFE, MATPRO)
%
% PURPOSE:
%   Computes the consistent linear stiffness matrix `K` at the undeformed
%   configuration (typically for small strains). This is used in:
%     - Dynamic simulations for Rayleigh damping (via β·K)
%     - ROM precomputation
%     - Simulations using linearized internal force evaluations
%
%   If requested via the `OTHER_output.ComputeLinearStiffnessMatrixForDamping_beta`
%   flag, it also contributes to the damping matrix `Ddamp` using the computed stiffness.
%
% CONDITIONS FOR EXECUTION:
%   The stiffness matrix is computed only if one of the following is true:
%     - Small strain formulation is used and `DATA.PRECOMPUTE_ELASTIC_STIFFNESS_MATRIX == 1`
%     - `PARAMETERS.OnlyComputeMassAndStiffnessMatricesSmallStrains == 1`
%     - `PARAMETERS.INTERNAL_FORCES_USING_precomputed_Kstiff == 1`
%     - `OTHER_output.ComputeLinearStiffnessMatrixForDamping_beta ~= 0`
%
% INPUT:
%   DATA         : Main problem data structure, including mesh, flags, etc.
%   PARAMETERS   : Additional flags and parameters for execution context
%   OTHER_output : Struct to store output matrices (e.g., stiffness K)
%   OPERFE       : Finite element operators and matrices (mass, damping, etc.)
%   MATPRO       : Material properties (stiffness tensor, etc.)
%
% OUTPUT:
%   OTHER_output.K     : Computed consistent linear stiffness matrix
%   OPERFE.Ddamp       : (If requested) Damping matrix with stiffness contribution
%   OPERFE.KinternalFORCES_given : (Optional) Precomputed stiffness used in force evaluation
%
% REMARKS:
%   - The stiffness is computed using `KstiffSmallStrains`, relying on
%     `celastST` and `FgradST` obtained from zero displacements.
%   - If only mass and stiffness matrices are requested, the function returns early.
%   - Optional commented-out block supports precomputing `celastINI_Bst` for nonlinear stress splitting (CECM).
%
% AUTHOR:
%   J.A. Hernández Ortega, Feb-2025, Balmes 185, Barcelona
%   Comments by ChatGPT4, 13-May-2025
% REFERENCE:
%   See also: 
%     - StressesFromDisplacementsVAR
%     - KstiffSmallStrains
%     - /FE_HROM/LARGE_STRAINS/InputDataFunctions/PreProcessInputDataDyn1.m
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp3.mat')
end
OTHER_output = DefaultField(OTHER_output,'ComputeLinearStiffnessMatrixForDamping_beta',0) ; 
  
    if (DATA.SMALL_STRAIN_KINEMATICS == 1 && DATA.PRECOMPUTE_ELASTIC_STIFFNESS_MATRIX  == 1) || ...
            PARAMETERS.OnlyComputeMassAndStiffnessMatricesSmallStrains == 1 || PARAMETERS.INTERNAL_FORCES_USING_precomputed_Kstiff == 1 ||...
        OTHER_output.ComputeLinearStiffnessMatrixForDamping_beta ~=0
            
        % Determine Stiffness Matrix, small strains
        %d  = zeros(DATA.MESH.ndof,1) ;
        VAR.DISP =zeros(DATA.MESH.ndof,1) ;VARint_n =[] ;
        
        
        
        [~,celastST,FgradST,~] = StressesFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
        OTHER_output.K = KstiffSmallStrains(OPERFE,FgradST,DATA.MESH.ndim,celastST) ;
        
        if OTHER_output.ComputeLinearStiffnessMatrixForDamping_beta >0 
            % Damping matrix, contribution stiffness matrix, see 
            % % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/InputDataFunctions/PreProcessInputDataDyn1.m

            betaD = OTHER_output.ComputeLinearStiffnessMatrixForDamping_beta; 
            OPERFE.Ddamp =  OPERFE.Ddamp + betaD*OTHER_output.K  ; 
        end
        
        if PARAMETERS.OnlyComputeMassAndStiffnessMatricesSmallStrains == 1
             return
        end
        if PARAMETERS.INTERNAL_FORCES_USING_precomputed_Kstiff == 1
            OPERFE.KinternalFORCES_given = OTHER_output.K ;
        end
    end
    
% DATA = DefaultField(DATA,'CECM_ONLY_FOR_NONLINEAR_STRESSES',0) ;
% 
% if DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 1
%     VAR.DISP =zeros(DATA.MESH.ndof,1) ;VARint_n =[] ;
%     [~,celastST_ini,~,~] = StressesFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
%     
%     celastST_iniD = ConvertBlockDiag(celastST_ini) ; % Diagonal block matrix
%     OPERFE.celastINI_Bst= celastST_iniD*OPERFE.Bst;
%     
%     
% end

    
    
 %    