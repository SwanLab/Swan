function [VAR,celastST,FgradST,detFgrad] = StressesFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n)
%--------------------------------------------------------------------------
% [VAR, celastST, FgradST, detFgrad] = StressesFromDisplacementsVAR(OPERFE, VAR, MATPRO, DATA, VARint_n)
%
% PURPOSE:
%   Computes the Green-Lagrange strains and the associated 2nd Piola-Kirchhoff
%   (PK2) and 1st Piola-Kirchhoff (PK1) stresses from the displacement field
%   stored in VAR.DISP. This function is central to nonlinear finite element 
%   analysis, enabling the evaluation of internal stresses based on the 
%   current deformation state.
%
% INPUT:
%   OPERFE    - Finite element operators (strain-displacement matrix Bst, weights, etc.).
%   VAR       - Structure with current state variables (mainly VAR.DISP).
%   MATPRO    - Material property structure (contains constitutive model data).
%   DATA      - Structure with global analysis settings (e.g., small/large strain flag,
%               use of deformation gradient, internal variables toggle, etc.).
%   VARint_n  - Structure with internal variables at the previous time step (plasticity models).
%
% OUTPUT:
%   VAR       - Updated structure including:
%                  VAR.GLSTRAINS     : Green-Lagrange strains at Gauss points
%                  VAR.PK2STRESS     : 2nd Piola–Kirchhoff stresses
%                  VAR.PK1STRESS     : 1st Piola–Kirchhoff stresses
%   celastST  - Consistent elastic tangent moduli at Gauss points (Voigt notation)
%   FgradST   - Deformation gradient at Gauss points (vectorized form)
%   detFgrad  - Determinant of the deformation gradient (only for large strains)
%
% NOTES:
%   - Automatically switches between small and large strain kinematics.
%   - Compatible with incremental and total-form constitutive updates.
%   - Uses `PK2stress_Constitutive_ModelVAR` to evaluate stresses and tangent moduli.
%   - When DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 1, also updates the nonlinear
%     incremental stress variable via `NonLinearStress_Incre`.
%--------------------------------------------------------------------------

% Copy of StressesFromDisplacements.m
% [VAR,celastST,FgradST,detFgrad] = StressesFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA) ;
%[GLSTRAINS,PK2STRESS,celastST,FgradST,PoneST,detFgrad] = StressesFromDisplacements(OPERFE,d,MATPRO,DATA) ;

% See
if nargin == 0
    load('tmp3.mat')
end


ndof = size(OPERFE.Bst,2) ;
FgradST =  OPERFE.Bst*VAR.DISP(1:ndof,:)  ; %+ OPERFE.IDENTITY_F ;
if DATA.SMALL_STRAIN_KINEMATICS ==0
    for idim = 1:DATA.MESH.ndim
        LOCROWS = idim:DATA.MESH.ndim^2:size(FgradST,1) ;
        FgradST(LOCROWS,:) = 1+FgradST(LOCROWS,:) ;
    end
    % 3. Green-Lagrante strains at all Gauss points
    VAR.GLSTRAINS = StrainGreenLagrange(FgradST,DATA.MESH.ndim) ;
else
    if DATA.NO_USE_Deformation_gradient_in_Small_Strains ==0
        % In small strain kinematics, FgradST == gradU (because we haven't added the identity matrix )
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
    [VAR,celastST,detFgrad ]= PK2stress_Constitutive_ModelVAR(VAR,MATPRO,DATA,FgradST,VARint_n,OPERFE) ;
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
        [VAR,celastST,detFgrad ]= PK2stress_Constitutive_ModelVAR(VAR,MATPRO,DATA,FgradSTloc,VARint_n) ;
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

%  For the case in which DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 1
VAR = NonLinearStress_Incre(DATA,OPERFE,VAR,FgradST) ;


