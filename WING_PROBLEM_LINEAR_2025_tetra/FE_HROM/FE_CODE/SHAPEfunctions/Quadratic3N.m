function  [Ne BeXi] = Quadratic3N(xi) ; 
% Shape functions and derivatives for 3-node quadratic element 
    % https://www.gidhome.com/documents/referencemanual/PREPROCESSING/Mesh%20Menu/Element%20type

if nargin == 0 
    % Check that 
    xi = -1 ; 
 %   Ne = [1 0 0]; 
   % xi = 0 ; 
  %  Ne = [0 1 0] ; 
%     xi = 1 ; 
%     Ne = [0 0 1] ; 
end
 
% Matrix of shape functions
Ne =[-0.5*xi*(1-xi),  0.5*xi*(1+xi),  (1-xi)*(1+xi)     ]; 
% Matrix of the gradient of shape functions 
BeXi = [ xi - 1/2,xi + 1/2, -2*xi ] ; 
