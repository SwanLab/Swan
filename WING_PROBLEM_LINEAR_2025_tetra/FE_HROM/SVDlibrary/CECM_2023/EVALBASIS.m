function  [PHIk_y,dPHIk_y,POLYINFO]=     EVALBASIS(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO)
%--------------------------------------------------------------------------
% function [PHIk_y, dPHIk_y, POLYINFO] = EVALBASIS(xNEW, DATA, VAR_SMOOTH_FE, POLYINFO)
%
% PURPOSE:
%   Evaluates the values and (optionally) the gradients of the empirical 
%   basis functions at the cubature points `xNEW`. This function supports 
%   two evaluation strategies:
%
%   1. **Direct Polynomial Fit (DATA.Method_Evaluation_Basis_Integrand == 1)**:
%      - Computes basis function values using fitted polynomial coefficients 
%        associated with the element containing each point.
%
%   2. **Analytical Evaluation (DATA.Method_Evaluation_Basis_Integrand == 2)**:
%      - Computes basis and gradient values analytically, based on shape functions.
%
%   The function is a key component in the **CECM (Continuous Empirical Cubature Method)** 
%   since it allows interpolation and evaluation of the reduced integrand at
%   arbitrary positions of integration points, which are updated during optimization.
%
% INPUT:
%   - xNEW           : [n x ndim] Coordinates of evaluation points (CECM points).
%   - DATA           : Struct with options and integrand description:
%                      .Method_Evaluation_Basis_Integrand : 1 = polynomial fit,
%                                                           2 = analytical.
%                      .Integrand : (if method 2) struct with gradient flag, etc.
%   - VAR_SMOOTH_FE  : Struct with mesh, polynomial orders, connectivities, etc.
%   - POLYINFO       : Struct with:
%                       .setElements : elements associated to xNEW,
%                       .TriangulationDelaunay : search structures for locating elements,
%                       .COEFFSpolynomial : polynomial coefficients (if available),
%                       .SCALING_VARIABLES : scaling and centering per element.
%
% OUTPUT:
%   - PHIk_y         : [r x n] Matrix of basis function values at xNEW.
%                      Each column corresponds to a reduced basis function φᵢ.
%   - dPHIk_y        : Gradient of basis functions (if available/applicable).
%   - POLYINFO       : Updated with new evaluation/triangulation info.
%
% DEPENDENCIES:
%   - EvaluateBasisFunctionDIRECTFIT2023
%   - EvaluateBasisFunctionANALYTICAL2023
%   - DefaultField (utility to fill missing struct fields)
%
% EXAMPLES:
%   Used within MAKE1ZERO_1STEP and SPARSIFglo to compute the residual 
%   Φᵗ·w - b and adjust integration point positions.
%
% REFERENCE:
%   J.A. Hernández et al., "Continuous Empirical Cubature Method for Data Compression
%   in Computational Mechanics", *Computer Methods in Applied Mechanics and Engineering*, 2024.
%
% AUTHOR:
%   Joaquín A. Hernández (UPC-CIMNE), 2023.
%--------------------------------------------------------------------------

% See comments in 
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/Poly1D/README_Quadrature1D.mlx
 
DATA = DefaultField(DATA,'Method_Evaluation_Basis_Integrand',1) ;      
POLYINFO = DefaultField(POLYINFO,'setElements',[]) ; 

POLYINFO = DefaultField(POLYINFO,'TriangulationDelaunay',cell(size(xNEW,1),1)) ; 

if DATA.Method_Evaluation_Basis_Integrand == 1 
    % Polynomial fitting (at each finite element)
    % --------------------------------------------
   % [PHIk_y,dPHIk_y,POLYINFO]=  EvaluateBasisFunctionDIRECTFIT(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO) ; 
    
    [PHIk_y,dPHIk_y,POLYINFO]=  EvaluateBasisFunctionDIRECTFIT2023(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO) ; 
    
elseif DATA.Method_Evaluation_Basis_Integrand == 2 
    DATA.Integrand.EVALUATE_GRADIENT = 1; 
     [PHIk_y,dPHIk_y,POLYINFO]=  EvaluateBasisFunctionANALYTICAL2023(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO) ; 
end
