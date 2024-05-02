function [N,Nxi,Neta] = ShapeFunc(pospg) 
% [N,Nxi,Neta] = ShapeFunc(pospg) 
% N, Nxi, Neta: matrices storing the values of the shape functions on the Gauss points
%               of the reference element
%               Each row concerns to a Gauss point
% pospg:        coordinates of Gauss points in the reference element
%

xi = pospg(:,1); eta = pospg(:,2); 
N    = [(1-xi).*(1-eta)/4, (1+xi).*(1-eta)/4, (1+xi).*(1+eta)/4, (1-xi).*(1+eta)/4]; 
Nxi  = [(eta-1)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4]; 
Neta = [(xi-1)/4, -(1+xi)/4,   (1+xi)/4,  (1-xi)/4 ]; 
