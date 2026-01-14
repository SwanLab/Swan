function  Fdis_e = FdisElem(tracBe,weig,shapef,dershapef,Xe) 
% Given 
% tracBe: Nodal values of the prescribed tractions at boundary element e (nnodeEb x1) (for a given spatial direction "i")
% weig :   Vector of Gauss weights (1xngaus), % shapef:   Array with the   shape functions at each Gauss point (ngaus x nnodeEb ), % dershapef:   Array with the derivatives of shape functions, with respect to
% element coordinates (ndimB x nnodeEb x ngaus), % Xe: Global coordinates of the nodes of the element,  
% % this function returns the element boundary tractions vector due to distributed loads  (Fdis_e)
if nargin == 0
    load('tmp.mat')
end
ngaus = length(weig) ; nnodeEb = size(Xe,2)  ; 
ndimB = size(dershapef,1) ;  Fdis_e = zeros(nnodeEb,1) ; 
for  g = 1:ngaus
    % Matrix of derivatives for Gauss point "g"
    BeXi = dershapef(:,:,g) ; 
    % Matrix of shape functions at point "g"
    Ne = shapef(g,:) ;    
    %%% Change of coordinates  
    Se = ChangeCoordBnd(Xe,ndimB) ; 
    %%% 
    % Jacobian Matrix 
    Je = Se*BeXi' ;  
    % JAcobian 
    detJe = det(Je) ;    
    %
    Fdis_e = Fdis_e + weig(g)*detJe*(Ne'*Ne)*tracBe ; 
end