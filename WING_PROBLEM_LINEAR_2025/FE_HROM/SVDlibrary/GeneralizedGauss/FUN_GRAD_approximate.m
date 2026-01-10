function [g,Grad_g] = FUN_GRAD_approximate(xNEW,DATAFITTING)

ndim = size(xNEW,2) ; 

if ndim == 1
    [g,Grad_g] = FUN_GRAD_approximate1D(xNEW,DATAFITTING) ; 
elseif ndim == 2 
    [g,Grad_g] = FUN_GRAD_approximate2D(xNEW,DATAFITTING) ; 
elseif ndim == 3 
    [g,Grad_g] = FUN_GRAD_approximate3D(xNEW,DATAFITTING) ; 
else
    error('Option not implemented')
end
