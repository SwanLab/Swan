function [VAR,celastST,FgradST,detFgrad] = StressesFromDisplacementsVARincre(OPERFE,VAR,MATPRO,DATA,VARint_n)
% ---------------------------------------------------------------------------------------------------
% FUNCTION: StressesFromDisplacementsVARincre
% ---------------------------------------------------------------------------------------------------
% PURPOSE:
%   Computes the **incremental nonlinear stress response** (second and first Piola–Kirchhoff tensors)
%   given a displacement field `VAR.DISP`. It supports both large strain and small strain kinematics,
%   and includes optional treatment of internal variables (for future plasticity extensions).
%
%   This function is a variation of `StressesFromDisplacementsVAR`, adapted specifically to:
%     - Handle **nonlinear strain–stress relationships**
%     - Support **incremental constitutive models** through `PK2stress_Constitutive_ModelVARincre`
%     - Allow postprocessing of multiple snapshots if needed
%
% USAGE:
%   [VAR, celastST, FgradST, detFgrad] = ...
%       StressesFromDisplacementsVARincre(OPERFE, VAR, MATPRO, DATA, VARint_n)
%
% INPUTS:
%   - OPERFE     : Structure with FE operators (e.g., strain–displacement matrices `Bst`).
%   - VAR        : Structure containing current displacement field and (to be filled) stress/strain data.
%   - MATPRO     : Structure with material properties for the constitutive law.
%   - DATA       : Structure with simulation flags and options (strain regime, Jacobian computation, etc.).
%   - VARint_n   : (Optional) Structure with internal variables from the previous time step.
%
% OUTPUTS:
%   - VAR        : Updated variable structure including:
%                   * `GLSTRAINS`  – Green-Lagrange strains or small strain tensor
%                   * `PK2STRESS`  – Second Piola–Kirchhoff stress tensor
%                   * `PK1STRESS`  – First Piola–Kirchhoff stress tensor
%   - celastST   : Consistent material tangent matrix (if requested)
%   - FgradST    : Deformation gradient at Gauss points (only in large strain case)
%   - detFgrad   : Determinant of deformation gradient (if needed by the constitutive model)
%
% FUNCTIONALITY:
%   - Computes the deformation gradient `F` from `Bst * u`.
%   - For large strain: F := I + ∇u; Green-Lagrange strains E = 0.5(FᵗF - I)
%   - For small strain: ∇u directly used as ε; identity may be skipped or enforced.
%   - Depending on `DATA.SMALL_STRAIN_KINEMATICS`, uses:
%       * `StrainGreenLagrange`
%       * `StrainGreenLagrange_small`
%   - Constitutive model is called via `PK2stress_Constitutive_ModelVARincre`
%   - Supports snapshot mode (loop over multiple displacement vectors), with internal variable update.
%
% NOTES:
%   - Internal variable logic is prepared but not implemented for large strain plasticity.
%   - If `FgradST` is empty (i.e., small strain regime), PK1 stresses are set to `[]`.
%
% REFERENCES:
%   - Theory and usage in:
%     /TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
%
% AUTHOR:
%   Joaquín A. Hernández, UPC – CIMNE  
%   Date: 6-Jan-2024, Balmes 185, Barcelona
%   Commented by ChatGPT.-4, 12-May-2025
%
% DEPENDENCIES:
%   - StrainGreenLagrange
%   - StrainGreenLagrange_small
%   - PK2stress_Constitutive_ModelVARincre
%   - PK1stress
%   - DefaultField
%
% ---------------------------------------------------------------------------------------------------

% Copy of StressesFromDisplacementsVAR.m 
% Adapted to determine "incremental" stresses (Nonlinear)
%  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
% 6-Jan-2024, Balmes, Barcelona 
%--------------------------------------
if nargin == 0
    load('tmp3.mat')
end

ndof = size(OPERFE.Bst,2) ;
FgradST =  OPERFE.Bst*VAR.DISP(1:ndof,:)  ; %+ OPERFE.IDENTITY_F ;
if DATA.SMALL_STRAIN_KINEMATICS ==0
  %  error('Option not implemented')
    for idim = 1:DATA.MESH.ndim
        LOCROWS = idim:DATA.MESH.ndim^2:size(FgradST,1) ;
        FgradST(LOCROWS,:) = 1+FgradST(LOCROWS,:) ;
    end
    % 3. Green-Lagrante strains at all Gauss points
    VAR.GLSTRAINS = StrainGreenLagrange(FgradST,DATA.MESH.ndim) ;
else
    if DATA.NO_USE_Deformation_gradient_in_Small_Strains ==0
        % In small strain kinematics, FgradST == gradU (because we haven't added the identity matrix )
       % error('Option not implemented')
        [VAR.GLSTRAINS,FgradST] = StrainGreenLagrange_small(FgradST,DATA.MESH.ndim) ;
    else
        % 8-Feb-2022. See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
        VAR.GLSTRAINS = FgradST ; % Small strains are used, from the beginning
        FgradST = [] ;
    end
end
% 4. 2nd Piola-Kirchhoff stresses at all Gauss Points
DATA.CALC_CTANG = 1 ;
if size(VAR.GLSTRAINS,2) == 1
    %  [PK2STRESS,celastST,detFgrad ]= PK2stress_Constitutive_Model(GLSTRAINS,MATPRO,DATA,FgradST) ;
    [VAR,celastST,detFgrad ]= PK2stress_Constitutive_ModelVARincre(VAR,MATPRO,DATA,FgradST,VARint_n) ;
else
    % This is used only for post-process purposes (not in the actual Newton-Raphson algorithm)
    % error('ADapt this option to the new framework (9-Feb-2022)')
    GLSTRAINS  = VAR.GLSTRAINS ;
    PK2STRESS = zeros(size(GLSTRAINS)) ;   celastST = [] ; detFgrad = [] ;
    if ~isempty(VARint_n)
        FFF = fieldnames(VARint_n)   ;
     %   VARINTtime = [] ;
        for ivar = 1:length(FFF)
            FLOC = FFF{ivar} ;
            VARint_n.(FLOC) = VAR.(FLOC) ;
            % VARINTtime.(FLOC) = zeros(size(VAR.(FLOC),1),size(GLSTRAINS,2)) ;
        end
    end
    % VAR.GLSTRAINS
    DATA.CALC_CTANG = 0 ;
    DATA.kiter = 2;
    for  itime = 1:size(GLSTRAINS,2)
        VAR.GLSTRAINS = GLSTRAINS(:,itime) ;
        if isempty(FgradST)
            FgradSTloc = [] ;
        else
            FgradSTloc = FgradST(:,itime) ;
        end
        [VAR,celastST,detFgrad ]= PK2stress_Constitutive_ModelVARincre(VAR,MATPRO,DATA,FgradSTloc,VARint_n) ;
        PK2STRESS(:,itime) = VAR.PK2STRESS ;
        
        % Internal variables
        if ~isempty(VARint_n)
            FFF = fieldnames(VARint_n)   ;
            for ivar = 1:length(FFF)
                FLOC = FFF{ivar} ;
                VARint_n.(FLOC) = VAR.(FLOC) ;
                
            end
        end
        
        
    end
    VAR.PK2STRESS = PK2STRESS ;
end

if isempty(VAR.PK2STRESS)
    VAR.PK1STRESS = [] ;
else
    % 5. 1st Piola-Kirchhoff stresses at all Gauss Points
    if ~isempty(FgradST)
        VAR.PK1STRESS = PK1stress(VAR.PK2STRESS,FgradST,DATA.MESH.ndim) ;
    else
        % This means that we are in the small strain regime.
        VAR.PK1STRESS = [] ;% PK2STRESS ;
    end
    
end


