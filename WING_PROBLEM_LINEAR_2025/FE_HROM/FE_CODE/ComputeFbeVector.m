function Fbe = ComputeFbeVector(fe,weig,shapef,dershapef,Xe) ; 
% Given 
% fe: Nodal values of the body force    (nnodeE*ndim x1)
% weig :   Vector of Gauss weights (1xngaus)
% shapef:   Array with the   shape functions at each Gauss point (ngaus x nnodeE )
% dershapef:   Array with the derivatives of shape functions, with respect to
% element coordinates (ndim x nnodeE x ngaus)
% Xe: Global coordinates of the nodes of the element,  
% % this function returns the element body force vector  Fbe
ndim = size(Xe,1) ; ngaus = length(weig) ; nnodeE = size(Xe,2)  ; 
Fbe = zeros(ndim*nnodeE,1) ; 
for  g = 1:ngaus
    % Matrix of derivatives for Gauss point "g"
    BeXi = dershapef(:,:,g) ; 
    % Matrix of shape functions at point "g"
    NeSCL = shapef(g,:) ; 
    % Matrix of shape functions at point "g" (for vector-valued fields)
    Ne = StransfN(NeSCL,ndim) ; 
    % Jacobian Matrix 
    Je = Xe*BeXi' ; 
    % JAcobian 
    detJe = det(Je) ;    
    %
    Fbe = Fbe + weig(g)*detJe*Ne'*(Ne*fe) ; 
end