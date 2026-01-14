function FbndE = ComputeFbndVector(qFLUXe,weig,shapef,dershapef,Xe) ; 
% Given 
% qFLUXe: Nodal values of the prescribed flux at boundary element e (nnodeEb x1)
% weig :   Vector of Gauss weights (1xngaus)
% shapef:   Array with the   shape functions at each Gauss point (ngaus x nnodeEb )
% dershapef:   Array with the derivatives of shape functions, with respect to
% element coordinates (ndimB x nnodeEb x ngaus)
% Xe: Global coordinates of the nodes of the element,  
% 
% this function returns the element boundary flux vector  FbndE
if nargin == 0
    load('tmp3.mat')
end

ngaus = length(weig) ; nnodeEb = size(Xe,2)  ; 
ndimB = size(dershapef,1) ; 
FbndE = zeros(nnodeEb,1) ; 

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
    FbndE = FbndE + weig(g)*detJe*(Ne'*Ne)*qFLUXe ; 
end
