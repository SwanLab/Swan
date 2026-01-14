function  [Ne BeXi] = Tetrahedra4N(xiV) ; 
% Shape functions and derivatives for 4-node tetrahedral 
xi = xiV(1) ; eta = xiV(2) ; zeta = xiV(3) ; 
% Matrix of shape functions
Ne = [1-eta-xi-zeta,eta, xi, zeta];

% Matrix of the gradient of shape functions 
BeXi = [-1,1, 0, 0 ;
        -1,0, 1, 0;
        -1,0, 0, 1];
    