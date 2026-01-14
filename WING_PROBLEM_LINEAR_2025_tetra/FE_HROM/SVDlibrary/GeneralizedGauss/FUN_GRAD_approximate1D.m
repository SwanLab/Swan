function [g,Grad_g] = FUN_GRAD_approximate1D(xNEW,DATAFITTING)

if nargin == 0
    load('tmp.mat')
end
nfun = length(DATAFITTING.spline_G) ; 
g= zeros(nfun,length(xNEW)); 
dg = zeros(size(g))  ; 
for i=1:nfun
    g(i,:) = fnval(DATAFITTING.spline_G{i},xNEW) ;
    dg(i,:) = fnval(DATAFITTING.spline_derG{i},xNEW) ;
end
Grad_g = {dg} ; 