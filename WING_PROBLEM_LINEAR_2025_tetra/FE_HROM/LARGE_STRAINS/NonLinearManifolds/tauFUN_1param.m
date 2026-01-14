function [tau,tauDER1,tauDER2] = tauFUN_1param(qL,DATA)
if nargin == 0
    load('tmp1.mat')
end
% THIS IS A FUNCTION TO EVALUATE THE MODAL COEFFICIENTS OF A REDUCED ORDER
% MODEL IN WHICH THE DISPLACEMENTS IS EVALUATED AS
% d_L = Phi*tau(qMASTER) = [PhiMASTER,PhiSLAVE]*[qSLAVE,f(qSLAVE)]
% Here qL = qMASTER is the sole generalized coordinate
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/NonLinearManifolds/BsplinesLeastSquares_fast.m
% (this is the function where the necessary inputs are defined )

% In turn, the nonlinear function f is determined by
% f = UU*diag(SS)*g(qL)

% DATA.sp = sp ;  % Function
% DATA.sp1 = sp1 ; % Derivatives function
% DATA.sp2 = sp2 ; % Derivatives function
% DATA.xmin = xmin ;
% DATA.xmax = xmax ;
% DATA.UleftSingular = UU ;
% DATA.UleftSingular = SS ;

% JAHO, 10th August 2025, Molinos Marfagones, Cartagena

f = zeros(length(DATA.SSingular),length(qL)) ;
df = f ; 
d2f = f; 
% Evaluate function
for islave = 1:length(DATA.SSingular)
    % We first evaluate the value and derivatives of the islave function
    [f(islave,:),df(islave,:),d2f(islave,:) ] = evaluate_spline_with_extrapolationFUN(DATA.sp{islave}, DATA.sp1{islave}, DATA.sp2{islave}, DATA.xmin,DATA.xmax,qL) ; 
end
% Thus 
tau = [qL;DATA.UleftSingular*(DATA.SSingular.*f)] ; 
tauDER1 = [ones(size(qL));DATA.UleftSingular*(DATA.SSingular.*df)] ; 
tauDER2 = [zeros(size(qL));DATA.UleftSingular*(DATA.SSingular.*d2f)] ; 



