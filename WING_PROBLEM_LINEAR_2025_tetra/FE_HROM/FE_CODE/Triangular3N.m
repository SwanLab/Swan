function  [Ne BeXi] = Triangular3N(xiV) ; 
% Shape functions and derivatives for triangles
xi = xiV(1) ; eta = xiV(2) ; 
% Matrix of shape functions
Ne =[1-xi-eta, xi, eta ]; 
% Matrix of the gradient of shape functions 
BeXi = [-1 1  0
        -1 0  1] ; 
