function  [Ne BeXi] = Linear2N(xi) ; 
% Shape functions and derivatives for 2-node linear element 
 
% Matrix of shape functions
Ne =0.5*[(1-xi)  (1+xi)  ]; 
% Matrix of the gradient of shape functions 
BeXi = 0.5*[-1 1] ; 
