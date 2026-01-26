function [g,Grad_g] = LocEvalBasFunAnalytical2023(xLOC,DATA,VAR_SMOOTH_FE,POLYINFO)
% Evaluation of basis functions and their gradients at  a given point xLOC
% See EvaluateBasisFunctionANALYTICAL.m
if nargin == 0
    load('tmp.mat')
end
% Evaluation of integrand functions and their gradients at xLOC
[f,df] = feval(DATA.Integrand.NameFunctionGenerate ,xLOC,DATA.Integrand) ;
if ~iscell(df)
    df = {df} ;
end
g =(f*VAR_SMOOTH_FE.VSinv) ;
Grad_g = cell(size(df)) ;

if DATA.Integrand.EVALUATE_GRADIENT == 1
    for idim = 1:length(df)
        Grad_g{idim} = df{idim}*VAR_SMOOTH_FE.VSinv ;
    end
end

 