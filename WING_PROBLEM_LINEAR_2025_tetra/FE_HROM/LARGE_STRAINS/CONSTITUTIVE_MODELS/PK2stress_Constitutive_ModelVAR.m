function [VAR,celasST,detFgrad ]= PK2stress_Constitutive_ModelVAR(VAR,MATPRO,DATA,FgradST,VARint_n,OPERFE)
%--------------------------------------------------------------------------
% PK2stress_Constitutive_ModelVAR
% --------------------------------
% Copy of  PK2stress_Constitutive_Model.m
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
% JAHO, 28-Oct-2024, UPC, Terrassa
%
% PURPOSE
%   Evaluate, at Gauss points, the constitutive response and (elastic)
%   consistent tangent for the current strain state, returning:
%     • PK2 stresses (VAR.PK2STRESS),
%     • elastic tangent (celasST) in Voigt form,
%     • (optionally) detF = det(F) when FgradST is supplied.
%   The routine supports small-strain elasticity, Neo-Hookean hyperelasticity,
%   small-strain J2 plasticity with history update, and a small-strain
%   isotropic damage law with linear hardening.
%
% INPUTS
%   VAR       : State with fields used here:
%               • GLSTRAINS  – strain measure at Gauss points (definition
%                 consistent with DATA.TYPE_CONSTITUTIVE_MODEL_ALL)
%               • (updated inside on plastic/damage paths)
%   MATPRO    : Material parameters (can be per-Gauss-point arrays/structs).
%   DATA      : Run/model controls, including:
%               • MESH.ndim, MESH.ngaus (optional, for reporting)
%               • PROCEED_WITH_NEGATIVE_JACOBIANS (0/1)
%               • TYPE_CONSTITUTIVE_MODEL_ALL ∈ {
%                     'SMALL_STRAINS_ELASTIC',
%                     'NeoHookean',
%                     'SMALL_STRAINS_J2_PLASTICITY',
%                     'SMALL_STRAINS_ISOTROPIC_DAMAGE_MODEL_LINEAR_HARD' }
%   FgradST   : (optional) deformation gradients at Gauss points. Required for
%               large-strain hyperelastic path ('NeoHookean').
%   VARint_n  : (optional) internal variables at previous step (plasticity/
%               damage models will read & update).
%   OPERFE    : (optional) not used here; reserved for future extensions.
%
% OUTPUTS
%   VAR       : Augmented/updated with:
%               • PK2STRESS  – second Piola–Kirchhoff stress at Gauss points
%               • (plastic/damage) updated internal variables
%               • NEGATIVE_JACOBIANS – boolean flag if detF≤0 detected
%   celasST   : Elastic tangent operator at Gauss points (Voigt).
%   detFgrad  : det(F) at Gauss points (empty if FgradST is not provided).
%
% BEHAVIOR & SAFEGUARDS
%   • If FgradST provided, detFgrad = Determinant_Fgrad(FgradST, DATA.MESH.ndim).
%     - If any detF≤0 and DATA.PROCEED_WITH_NEGATIVE_JACOBIANS==0:
%         * Warn, print offending elements (when MESH.ngaus available),
%         * Return empty stress/tangent and exit early.
%     - If detF≤0 and proceeding is allowed, set VAR.NEGATIVE_JACOBIANS=true,
%       continue the constitutive evaluation (use with care).
%   • Model selection by DATA.TYPE_CONSTITUTIVE_MODEL_ALL:
%       'SMALL_STRAINS_ELASTIC'
%          → [VAR.PK2STRESS, celasST] = SmallStrainLargeRotations(GLSTRAINS, MATPRO)
%       'NeoHookean'
%          → [VAR.PK2STRESS, celasST] = NeoHookStressStrain(detFgrad, GLSTRAINS, MATPRO, DATA)
%       'SMALL_STRAINS_J2_PLASTICITY'
%          → [VAR, celasST]           = SmallStrainJ2Plasticity(VAR, MATPRO, DATA, VARint_n)
%       'SMALL_STRAINS_ISOTROPIC_DAMAGE_MODEL_LINEAR_HARD'
%          → [VAR, celasST]           = DamageSmallStrainLargeRotations(VAR, GLSTRAINS, MATPRO, VARint_n)
%     Unrecognized options raise an error.
%
% NOTES
%   • GLSTRAINS must match the kinematic assumption of the selected model.
%   • Plasticity/damage paths update VAR’s internal variables in-place.
%   • Hyperelastic Neo-Hookean path requires both FgradST and detFgrad.
%   • When early-exiting due to negative Jacobians (and prohibition flag),
%     upper layers should interpret empty outputs as “invalid state”.
%
% DEPENDENCIES
%   Determinant_Fgrad, large2smallREP,
%   SmallStrainLargeRotations, NeoHookStressStrain,
%   SmallStrainJ2Plasticity, DamageSmallStrainLargeRotations.
%
% AUTHOR
%   Joaquín A. Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%
% COMMENTS UPDATED BY: ChatGPT-5
%   Date of update: 07-NOV-2025, Barcelona, Spain
%--------------------------------------------------------------------------

% -----------------------------------------------------------------------------

if nargin == 0
    load('tmp.mat')
elseif nargin == 3
    FgradST = [] ; VARint_n = [] ; OPERFE = [] ;
elseif nargin == 4
    VARint_n = [] ; OPERFE = [] ;
elseif nargin == 5
    OPERFE = [] ;
end
detFgrad = [] ;

% Determinant Gradient FgradST
if ~isempty(FgradST)
    detFgrad = Determinant_Fgrad(FgradST,DATA.MESH.ndim) ;
    
    IND_NEG = find(detFgrad<=0) ;
else
    IND_NEG = [] ;
end


%IND_NEG = [] ;

%DATA.PROCEED_WITH_NEGATIVE_JACOBIANS = 1;  %15-May_2024

 VAR.NEGATIVE_JACOBIANS = false; 
if ~isempty(IND_NEG)  && DATA.PROCEED_WITH_NEGATIVE_JACOBIANS == 0
    warning('Negative Jacobian...')
    
    if isfield(DATA.MESH,'ngaus')
        IND_ELEM = large2smallREP(IND_NEG,DATA.MESH.ngaus) ;
        disp(['Elements with Negative Jacobians'])
        disp(IND_ELEM')
        %   clipboard('copy',num2str(IND_ELEM'))
    end
    
    StwoST = [] ; celasST = [] ; detFgrad = [] ;
    
    
else
    
    if  ~isempty(IND_NEG)  %15-May_2024
        %   disp('Negative Jacobian...')
        fprintf('Negative Jacobian detected.\n');
        if isfield(DATA.MESH,'ngaus')
            IND_ELEM = large2smallREP(IND_NEG,DATA.MESH.ngaus) ;
            %    disp(['Elements with Negative Jacobians'])
            fprintf('%d ', IND_ELEM);
            fprintf('\n');
            
        end
        %  disp('We proceed with the calculations though... (set DATA.PROCEED_WITH_NEGATIVE_JACOBIANS = 0 for proceeding otherwise)')
    VAR.NEGATIVE_JACOBIANS = true; 
   
    end
    
    switch  DATA.TYPE_CONSTITUTIVE_MODEL_ALL
        case 'SMALL_STRAINS_ELASTIC'
            [VAR.PK2STRESS,celasST]= SmallStrainLargeRotations(VAR.GLSTRAINS,MATPRO) ;
        case 'NeoHookean'
            [VAR.PK2STRESS,celasST]= NeoHookStressStrain(detFgrad,VAR.GLSTRAINS,MATPRO,DATA) ;
            
        case 'SMALL_STRAINS_J2_PLASTICITY'
            [VAR,celasST]= SmallStrainJ2Plasticity(VAR,MATPRO,DATA,VARint_n) ;
            
        case 'SMALL_STRAINS_ISOTROPIC_DAMAGE_MODEL_LINEAR_HARD'
             [VAR,celasST]= DamageSmallStrainLargeRotations(VAR,VAR.GLSTRAINS,MATPRO,VARint_n) ;
            
        otherwise
            error('OPtion not implemented yet')
    end
    
    
    
    
end
