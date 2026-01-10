function [dN_dx detJ] = QuadElem_NX(X,dN_dXi)

if nargin ==0
    X = [-1 -1;
        1 -1;
        1  1
        -1  1];
    [N dN_dXi] = QuadElem_N_derN  ;
end



ndime = size(X,2);

ngaus = size(dN_dXi,1)/ndime ;
detJ = zeros(1,ngaus) ;
dN_dx = zeros(size(dN_dXi)) ;
for g =1:ngaus
    % Jacobian matrix
    iniROW = (g-1)*ndime+1; finROW = g*ndime ;
    dN_dXi_g = dN_dXi(iniROW:finROW,:) ;
    J_g = dN_dXi_g*X ;
    % Determinant of Jacobian matrix
    detJ(g) = det(J_g) ;
    dN_dx_g = J_g\dN_dXi_g ;
    dN_dx(iniROW:finROW,:) = dN_dx_g;
end
