function [VAR,VARINTtime,celastST,FgradST,detFgrad] = StressesFromDisplacementsVAR_EIFE(OPERFE,VAR,MATPRO,DATA,VARint_n)
% Copy of StressesFromDisplacementsVAR.m
% EIFE method
% JHAO, 24-MArch-2024
% Comments by ChatGPT-4
% %StressesFromDisplacementsVAR_EIFE Compute stress tensors and internal variables from strain history in EIFE
%
%   This function computes the 2nd Piola-Kirchhoff stress tensor and internal variables 
%   for all Gauss points of an EIFE element, using the history of displacements projected 
%   via the B-matrix. It supports both small strain and finite strain kinematics, and 
%   handles internal variable tracking via a generic material interface.
%
%   INPUT:
%     OPERFE      - Structure with operators for computing strains (e.g., Bst matrix)
%     VAR         - Structure containing the displacement field history (VAR.DISP)
%     MATPRO      - Material properties at Gauss points, including constitutive model info
%     DATA        - Control parameters and solver options (e.g., strain regime flags)
%     VARint_n    - (Optional) Initial state of internal variables per Gauss point
%
%   OUTPUT:
%     VAR         - Structure augmented with:
%                     - GLSTRAINS: Green–Lagrange strains (or small strain tensor)
%                     - PK2STRESS: 2nd Piola–Kirchhoff stress
%                     - PK1STRESS: 1st Piola–Kirchhoff stress (if Fgrad known)
%     VARINTtime  - Time evolution of internal variables (struct of arrays)
%     celastST    - Consistent tangent operator (only used in Newton loops)
%     FgradST     - Deformation gradient at Gauss points and time steps
%     detFgrad    - Determinant of Fgrad at Gauss points and time steps
%
%   PROCEDURE:
%   -----------------------------------------------------------------------------------------
%   1. Compute gradient of displacements (∇u or Fgrad - I) from B-matrix and DOFs
%   2. Compute strain tensor:
%        - Green–Lagrange for finite strains
%        - Small strain tensor or directly ∇u for small strain regime
%   3. Evaluate stress and internal variable update via:
%        PK2stress_Constitutive_ModelVAR
%   4. If multiple time steps are present:
%        - Loop over time steps
%        - Track time evolution of stresses and internal variables
%   5. Optionally compute PK1 stress if Fgrad is available
%
%   SPECIAL CASES:
%     - If `DATA.SMALL_STRAIN_KINEMATICS = 1`, small strain formulation is used.
%     - If `DATA.NO_USE_Deformation_gradient_in_Small_Strains = 1`, gradient is identity.
%     - Supports both incremental and total strain–stress evaluation workflows.
%
%   REMARKS:
%     - This routine is central to the stress recovery step in EIFE, and interfaces
%       directly with user-defined or library constitutive models through `MATPRO`.
%     - Internally, `DATA.CALC_CTANG` determines whether the consistent tangent
%       operator is calculated (needed in Newton-Raphson).
%
%   AUTHOR:
%     Joaquín A. Hernández Ortega (JAHO), UPC BarcelonaTech, 24-March-2024
%     Adapted from `StressesFromDisplacementsVAR.m` for use in the EIFE framework
%
%   SEE ALSO:
%     PK2stress_Constitutive_ModelVAR, StrainGreenLagrange, PK1stress,
%     MaterialPropertiesIntegPoints, StrainGreenLagrange_small


% See
if nargin == 0
    load('tmp.mat')
end


ndof = size(OPERFE.Bst,2) ;
FgradST =  OPERFE.Bst*VAR.DISP(1:ndof,:)  ; %+ OPERFE.IDENTITY_F ;
if DATA.SMALL_STRAIN_KINEMATICS ==0
    for idim = 1:DATA.MESH.ndim
        LOCROWS = idim:DATA.MESH.ndim^2:length(FgradST) ;
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
    [VAR,celastST,detFgrad ]= PK2stress_Constitutive_ModelVAR(VAR,MATPRO,DATA,FgradST,VARint_n) ;
else
    % This is used only for post-process purposes (not in the actual Newton-Raphson algorithm)
    % error('ADapt this option to the new framework (9-Feb-2022)')
     GLSTRAINS  = VAR.GLSTRAINS ;
    PK2STRESS = zeros(size(GLSTRAINS)) ;   celastST = [] ; detFgrad = [] ;
    if ~isempty(VARint_n)
        FFF = fieldnames(VARint_n)   ;
        VARINTtime = [] ;
        for ivar = 1:length(FFF)
            FLOC = FFF{ivar} ;
      %      VARint_n.(FLOC) = VAR.(FLOC) ;
            VARINTtime.(FLOC) = zeros(size(VARint_n.(FLOC),1),size(GLSTRAINS,2)) ;
        end
    else
        VARINTtime = [] ; 
    end
    % VAR.GLSTRAINS
    DATA.CALC_CTANG = 0 ;
    DATA.kiter = 2;
    VARLOC = []; 
    for  itime = 1:size(GLSTRAINS,2)
        VARLOC.GLSTRAINS = GLSTRAINS(:,itime) ;
        if isempty(FgradST)
            FgradSTloc = [] ;
        else
            FgradSTloc = FgradST(:,itime) ;
        end
        [VARLOC,celastST,detFgrad ]= PK2stress_Constitutive_ModelVAR(VARLOC,MATPRO,DATA,FgradSTloc,VARint_n) ;
        PK2STRESS(:,itime) = VARLOC.PK2STRESS ;
        
        
        % Internal variables
        if ~isempty(VARint_n)
            FFF = fieldnames(VARint_n)   ;
            for ivar = 1:length(FFF)
                FLOC = FFF{ivar} ;
                VARint_n.(FLOC) = VAR.(FLOC) ;
                VARINTtime.(FLOC)(:,itime) =  VAR.(FLOC) ;
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


