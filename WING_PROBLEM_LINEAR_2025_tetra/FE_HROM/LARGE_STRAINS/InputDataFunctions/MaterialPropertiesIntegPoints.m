function    [MATPRO,MESH,DATA,INICOND] = MaterialPropertiesIntegPoints(MESH,DATA,PROPMAT) 
% %MaterialPropertiesIntegPoints Assign material properties and initialize internal variables at Gauss points
%
%   This function assigns material properties at Gauss integration points for the selected
%   constitutive model in the simulation. It also initializes internal variables required
%   by inelastic models (e.g., plasticity), and updates the `DATA` structure accordingly.
%
%   INPUT:
%     MESH     - Mesh structure of the EIFE domain, including Gauss points and element info
%     DATA     - Structure containing simulation settings, including:
%                  - DATA.TYPE_CONSTITUTIVE_MODEL_ALL: model selector
%                  - DATA.MESH.nstrain: number of strain components per Gauss point
%     PROPMAT  - Structure with material parameters (e.g., Young's modulus, yield stress)
%
%   OUTPUT:
%     MATPRO   - Structure with constitutive parameters assigned at Gauss points
%     MESH     - (Possibly updated) mesh structure
%     DATA     - Updated simulation data with initialized model-specific fields
%     INICOND  - Structure containing initial values of internal variables, such as:
%                  - YieldStress
%                  - PlasticStrains
%                  - InternalVarStrain
%
%   PROCEDURE:
%   ---------------------------------------------------------------------------------------
%   Depending on the selected constitutive model:
%
%     - 'SMALL_STRAINS_ELASTIC':
%         * Calls `SmallStrainElasticityPROP` to assign linear elastic properties.
%         * Constructs a global block-diagonal elasticity tensor for all Gauss points.
%
%     - 'NeoHookean':
%         * Calls `NeoHookPROP` for compressible Neo-Hookean hyperelasticity.
%         * No internal variables are needed.
%
%     - 'SMALL_STRAINS_J2_PLASTICITY':
%         * Calls `SmallStrainJ2PlasticityPROP` to define elastoplastic properties.
%         * Initializes internal variables: `YieldStress`, `PlasticStrains`, and `InternalVarStrain`.
%         * Registers them in `DATA.ListFieldInternalVariables`.
%
%   REMARKS:
%     - This function prepares the material and internal state for stress updates 
%       (via PK2stress_Constitutive_ModelVAR) and strain integration across time steps.
%     - It is typically called once at the beginning of the simulation or during domain setup.
%
%   AUTHOR:
%     Joaquín A. Hernández Ortega (JAHO), UPC BarcelonaTech, March 2024
%
%   SEE ALSO:
%     SmallStrainElasticityPROP, NeoHookPROP, SmallStrainJ2PlasticityPROP,
%     PK2stress_Constitutive_ModelVAR

DATA.ListFieldInternalVariables = [] ; % Use this variable to specify the list of internal variables 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
INICOND = [] ; 
switch DATA.TYPE_CONSTITUTIVE_MODEL_ALL
    case 'SMALL_STRAINS_ELASTIC'
        MESH.nstrain=  DATA.MESH.nstrain  ; 
        [MATPRO] = SmallStrainElasticityPROP(MESH,DATA.typePROBLEM,PROPMAT) ;
        disp('Global elasticity matrix    ...')
        Cglo = DefineElastMatGLO_nw(MATPRO.celasglo,DATA.MESH.ngaus) ; % ; ...
        %  Cglo = ConvertBlockDiag(Cglo) ;
        MATPRO.celasglo = [] ;
        MATPRO.celasglo = Cglo ;
        
    case 'NeoHookean'
        [MATPRO] = NeoHookPROP(MESH,DATA.typePROBLEM,PROPMAT,DATA) ;
    %    disp('Global elasticity matrix    ...')
        %     Cglo = DefineElastMatGLO_nw(MATPRO.celasglo,DATA.MESH.ngaus) ; % ; ...
        %  Cglo = ConvertBlockDiag(Cglo) ;
        %    MATPRO.celasglo = [] ;
        %   MATPRO.celasglo = Cglo ;
    case 'SMALL_STRAINS_J2_PLASTICITY'
        % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx
        [MATPRO,DATA] = SmallStrainJ2PlasticityPROP(MESH,DATA.typePROBLEM,PROPMAT,DATA) ;
        INICOND.YieldStress =  MATPRO.sigmay_0 ;   % Initial condition 
        DATA.ListFieldInternalVariables = {'YieldStress','PlasticStrains','InternalVarStrain'} ; 
        ngausT = size(INICOND.YieldStress,1) ; 
        INICOND.PlasticStrains = zeros(ngausT*DATA.MESH.nstrain,1) ; 
        INICOND.InternalVarStrain = zeros(size(INICOND.YieldStress)) ; 
        
    case 'SMALL_STRAINS_ISOTROPIC_DAMAGE_MODEL_LINEAR_HARD'
        
        [DATA,INICOND,MATPRO] = DamageMaterialPropIntPoints(DATA,MESH,PROPMAT) ; 
      
        
        
        
  
        
    otherwise
        error('Option not implemented')
end
