function Kelem = ElementKMatrices(nelem,nnodeE,ndim,BCB,ngaus)
% Given the matrix BCB containing the stiffness matrices at all gauss
% points, ElementKMatrices returns the matrix of element stiffness matrix,
% obtained by summing up the contribution of each Gauss point.
%dbstop('6')
if nargin==0
    BCB = [1 2 3 4 5 6 1 2 3 4 5 6]' ; nelem = 2 ; nnodeE = 3; ndim=1 ; ngaus=2;
    BCB = [BCB BCB BCB];
end

Kelem = zeros(nelem*nnodeE*ndim,nnodeE*ndim) ;
indREF = 1:nnodeE*ndim ;
ROWSref = repmat(indREF,nelem,1);
for g = 1:ngaus
    gini = (g-1)*nnodeE*ndim+1 ;
    indPOS = (gini:nnodeE*ndim*ngaus:nnodeE*ngaus*ndim*nelem)'-1 ;
    ROWS = bsxfun(@plus,ROWSref,indPOS);  ROWS =ROWS' ;
    ROWS = ROWS(:) ;
    Kelem = Kelem + BCB(ROWS,:) ;
end

end