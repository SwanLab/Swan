function   [STRESSES_FINE,INTERNAL_VAR,VMSTRESS] = ReconstStressesEIFE_exactBUBcorot(EIFE_prop,dClocINCREe_time,DATA,...
    lambdaLENe,PROPMAT ) 
% Reconstruction of stresses and internal variables by solving integrating
% the constitutive equations along time 
% Adaptation of ReconstStressesEIFE_exactBUB.m to the co-rotational
% approach 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
% JAHO, 6-Nov-2024, UPC, Terrassa. 
% Comments by CHATGPT-4
%ReconstStressesEIFE_exactBUBcorot Reconstruct fine-scale stresses using a co-rotational EIFE approach with bubble modes
%
%   This function computes the stress field and internal variables for a fine-scale
%   EIFEM domain using a history of local coarse-scale displacements, including
%   contributions from bubble modes. It adapts the standard reconstruction procedure
%   to account for large rotations via a co-rotational framework.
%
%   INPUT:
%     EIFE_prop         - Structure containing precomputed operators, mesh, and info for the EIFE domain
%     dClocINCREe_time  - Local coarse-scale incremental displacements (boundary + bubble DOFs) over time
%     DATA              - Simulation and postprocessing options (including material model settings)
%     lambdaLENe        - Scaling factor for the physical size of the domain
%     PROPMAT           - Material model database (with stress update routines)
%
%   OUTPUT:
%     STRESSES_FINE     - 2nd Piola–Kirchhoff stress tensor values at all Gauss points and time steps
%     INTERNAL_VAR      - Selected internal variable (e.g., strain energy, plastic strain) if defined
%     VMSTRESS          - Von Mises stress computed from the tensor at all Gauss points and time steps
%
%   MAIN STEPS:
%   ----------------------------------------------------------------------------------
%   1. Determine material properties at Gauss points via `MaterialPropertiesIntegPoints`
%   2. Build the strain projection operator:
%        B̂ = (1/λ)·(B_fine × Φ_DEF × P_down_DEF)
%   3. Feed the history of displacements into the constitutive model using:
%        `StressesFromDisplacementsVAR_EIFE`
%   4. Extract:
%      - STRESSES_FINE: full tensor history
%      - INTERNAL_VAR: selected internal variable, if requested
%      - VMSTRESS: Von Mises stress magnitude
%
%   REMARKS:
%   - The stresses are computed in the *reference configuration* (unrotated)
%     because the final transformation to global (rotated) configuration is not yet implemented.
%   - The function assumes small strains and uses PK2 stress output unless specified otherwise.
%   - It supports multiple internal variables, though only one can be selected for output.
%
%   TODO:
%   - Implement rotation of stress tensors to the global configuration.
%   - Generalize to support 2D plane strain rotation transformation.
%
%   AUTHOR:
%     Joaquín A. Hernández Ortega (JAHO), UPC BarcelonaTech, 6-Nov-2024
%     Adapted from ReconstStressesEIFE_exactBUB.m to support co-rotational EIFE
%
%   SEE ALSO:
%     StressesFromDisplacementsVAR_EIFE, MaterialPropertiesIntegPoints,
%     ReconstDisplacementsEIFEbubCOROT, VonMisesStressCOMP



%---------------------------------------------
if nargin == 0
    load('tmp1.mat')
end

% DETERMINATION OF VARIABLE "MATPRO" FOR ALL THE GAUSS POINTS OF THE EIF
% ELEMENT (in the parent configuration)
DATA.MESH.ngaus = size(EIFE_prop.MESH.posgp,2) ;
[MATPRO,MESH,DATA,INICOND] = MaterialPropertiesIntegPoints(EIFE_prop.MESH,DATA,PROPMAT) ;


if isempty(DATA.ListFieldInternalVariables)
    VARint_n = [] ;
else
    DATA = DefaultField(DATA,'NumberInternalVariableToPrint',[length(DATA.ListFieldInternalVariables)]);
    for ivar = 1:length(DATA.ListFieldInternalVariables)
        FLOC = DATA.ListFieldInternalVariables{ivar};
        VARint_n.(FLOC) = INICOND.(FLOC)  ;
        DATA.STORE.VAR.(FLOC) = 1;
    end
end
% STRESSES FROM HISTORY OF DISPLACEMENTS


%  Stresses can be computed using standard FE stress routines using as input $\dClocINCREe{e}(t)$, and as ``stacked'' global B-matrix: 
% \BstFClocE{e}  = \dfrac{1}{\lambdaLENe{e}} (\BstFlocE{e} \PhiDEFallE{e}) \PdownsDEFe{e}

% EIFEoper.RECONSTRUCTION.STRAINS.coeff = PdownsDEF ;
% EIFEoper.RECONSTRUCTION.STRAINS.BASIS = OPERFE.Bst*PhiDEF ;


OPERFE.Bst = EIFE_prop.RECONSTRUCTION.STRAINS.BASIS*EIFE_prop.RECONSTRUCTION.STRAINS.coeff/lambdaLENe ; 



VAR.DISP = dClocINCREe_time ;

[VAR,VARINTtime,~,FgradST,detFgrad] = StressesFromDisplacementsVAR_EIFE(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
%STRESSES_FINE_REF = VAR.PK2STRESS;

% Change introduce 13-May-2025, see  
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN.mlx
ndim = size(MESH.COOR,2) ;


if   DATA.SMALL_STRAIN_KINEMATICS ==0
    if isempty(detFgrad)
        detFgrad = Determinant_Fgrad(FgradST,DATA.MESH.ndim) ;
    end
    CauchyStress = CauchyStressFromPK1_GEN(VAR.PK1STRESS,FgradST,detFgrad,ndim) ;
    
else
    CauchyStress = VAR.PK1STRESS ; 
end
 
VMSTRESS = VonMisesStressCOMP(CauchyStress,ndim,DATA) ;

%VMSTRESS = VonMisesStressCOMP(STRESSES_FINE_REF,ndim,DATA) ;

if ~isempty(DATA.ListFieldInternalVariables)
    NameInternalVariable = DATA.ListFieldInternalVariables{DATA.NumberInternalVariableToPrint} ;
    INTERNAL_VAR = VARINTtime.(NameInternalVariable)  ;
else
    INTERNAL_VAR = [] ;
end


% if EIFE_prop.MESH.nstrain == 4
%     
%     STRESSES_FINE =zeros(size(STRESSES_FINE_REF)) ;
%     if length(ROTATIONmat) == 1
%         ROTATION_STRESSES = eye(3) ;
%     else
%         ROTATION_STRESSES =    RotateStress2Dplanestrain(ROTATIONmat(1,1),ROTATIONmat(2,1)) ;
%     end
%     
%     for itime = 1:size(STRESSES_FINE_REF,2)
%         STRESSES_FINE_loc =  reshape(STRESSES_FINE_REF(:,itime),EIFE_prop.MESH.nstrain,[])  ;
%         
%         STRESSES_FINE_loc(1:3,:) = ROTATION_STRESSES*STRESSES_FINE_loc(1:3,:) ;
%         
%         % GID CONVENTION FOR PRINTING
%         %   STRESSES_FINE = STRESSES_FINE([1 2 4 3],:) ;
%         STRESSES_FINE(:,itime) = STRESSES_FINE_loc(:);
%         %  ngaus =size(EIFE_prop.MESH.posgp,2) ;
%         %  nstrain = (EIFE_prop.MESH.nstrain) ;
%         
%         %  STRESSES_FINE_loc =  reshape(STRESSES_FINE_loc,nstrain*ngaus,[])  ;
%         
%     end
%     
%     
% else
    % ROTATION OF FINE-SCALE STRESSES NOT IMPLEMENTED YET
    STRESSES_FINE = VAR.PK2STRESS ;
    disp('Rotation of stresses not implemented yet ')
%end



