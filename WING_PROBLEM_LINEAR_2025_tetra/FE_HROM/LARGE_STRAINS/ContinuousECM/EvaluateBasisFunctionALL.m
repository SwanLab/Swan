function  [PHIk_y,dPHIk_y,POLYINFO]=     EvaluateBasisFunctionALL(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO)
% See comments in 
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/09_ASSESSMENTpaper/Poly1D/README_Quadrature1D.mlx
 
DATA = DefaultField(DATA,'Method_Evaluation_Basis_Integrand',1) ;      

if DATA.Method_Evaluation_Basis_Integrand == 1 
    % Polynomial fitting (at each finite element)
    % --------------------------------------------
    [PHIk_y,dPHIk_y,POLYINFO]=  EvaluateBasisFunctionDIRECTFIT(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO) ; 
elseif DATA.Method_Evaluation_Basis_Integrand == 2 
    DATA.Integrand.EVALUATE_GRADIENT = 1; 
     [PHIk_y,dPHIk_y,POLYINFO]=  EvaluateBasisFunctionANALYTICAL(xNEW,DATA,VAR_SMOOTH_FE) ; 
end
