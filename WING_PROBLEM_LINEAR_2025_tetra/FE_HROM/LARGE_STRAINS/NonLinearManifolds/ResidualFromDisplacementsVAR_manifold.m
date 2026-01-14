function  [VAR,celastST,FgradST,detFgrad] =  ResidualFromDisplacementsVAR_manifold(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
%--------------------------------------------------------------------------
% [VAR, celastST, FgradST, detFgrad] = ResidualFromDisplacementsVAR(OPERFE, VAR, MATPRO, DATA, VARint_n)
%
% PURPOSE:
%   Computes the internal forces and residual vector from the current 
%   displacement field using the constitutive model. This is the core routine 
%   to evaluate equilibrium equations in nonlinear solid mechanics problems.
%
% INPUT:
%   OPERFE    - Structure containing finite element operators (strain-displacement 
%               matrices, integration weights, etc.).
%   VAR       - Structure with state variables, especially VAR.DISP (nodal displacements).
%   MATPRO    - Material properties and constitutive law definitions.
%   DATA      - Structure with control flags and analysis parameters (e.g., nonlinear vs linear mode).
%   VARint_n  - Structure with internal variables at previous time step (plasticity, history-dependent models).
%
% OUTPUT:
%   VAR       - Updated structure with:
%                   VAR.PK1STRESS: First Piola-Kirchhoff stress tensor (vectorized)
%                   VAR.PK2STRESS: Second Piola-Kirchhoff stress tensor
%                   VAR.FINT: Internal force vector
%                   VAR.RESID: Residual vector (FINT - FEXT)
%   celastST  - Elastic tangent operator (Voigt notation) at all Gauss points
%   FgradST   - Deformation gradient (vectorized form) at all Gauss points
%   detFgrad  - Determinant of deformation gradient (used in nonlinear analysis)
%
% NOTES:
%   - Compatible with both small and large strain regimes.
%   - Supports precomputed stiffness matrix option for linear systems.
%   - Core function for Newton–Raphson iterations in nonlinear analysis.
%--------------------------------------------------------------------------

% This is a copy of ResidualFromDisplacements.m.
% [VAR,celastST,Fint,FgradST,detFgrad] =...
%        ResidualFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA) ;

%The inputs are not
% specified explicitily, nor the outputs.
% See  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/PLASTICITY_HROM/README_plasticity.mlx


if nargin == 0
    load('tmp3.mat')
end



% Stresses from displacements
%[GLSTRAINS,PK2STRESS,celastST,FgradST,PoneST,detFgrad] = StressesFromDisplacements(OPERFE,d,MATPRO,DATA) ;
if DATA.INTERNAL_FORCES_USING_precomputed_Kstiff ==0
    [VAR,celastST,FgradST,detFgrad] = StressesFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
else
    celastST = [] ; FgradST = [] ; detFgrad = [] ;
end


if isempty(VAR.PK2STRESS) && DATA.INTERNAL_FORCES_USING_precomputed_Kstiff ==0
    VAR.PoneST = [] ; Fint = [] ; VAR.RESID = [] ;
else
    OPERFE = DefaultField(OPERFE,'KinternalFORCES_given',[]) ;
    % 6.1. Internal forces
    if isempty(OPERFE.KinternalFORCES_given)
     %   if ~isfield(OPERFE,'BstW')           
        VAR.FINT = InternalForces(OPERFE,VAR.PK1STRESS,VAR.PK2STRESS,DATA,VAR) ;
      %  else 
            % 3-Jan-2024, see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/106_ECMtailoredTOmodes/01_VIABILITY.mlx
            % This option proved to be unreliable 
       %  VAR.FINT = InternalForces_TAILOREDWEIGHTS(OPERFE,VAR.PK1STRESS,VAR.PK2STRESS,DATA) ;
       % end
    else
        % FOR the dynamic mode decomposition, see
        %/home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/DynamicModeDecomposition/01_ForcedVibration_2D/RunDMD_aux.mlx
        
        % Internal forces are assumed to be linear with displacements (and given)
        if ~isstruct(OPERFE.KinternalFORCES_given)
            VAR.FINT = OPERFE.KinternalFORCES_given*VAR.DISP ;
        else
            iMATRIX = OPERFE.KinternalFORCES_given.INTERVALS(DATA.istep) ;
            Kloc = OPERFE.KinternalFORCES_given.MATRICES{iMATRIX} ;
            VAR.FINT = Kloc*VAR.DISP ;
        end
    end
    
    % 6.2. Residual
    VAR.RESID  = VAR.FINT- VAR.FEXT;
end


 