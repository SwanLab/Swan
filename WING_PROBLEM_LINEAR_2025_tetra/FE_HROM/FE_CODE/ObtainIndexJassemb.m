function j = ObtainIndexJassemb(B,ndim)
%dbstop('3')
if nargin == 0
    load('tmp1.mat')
end

B  = repmat(B,1,ndim );
nelem = length(B)/ndim ;
J = zeros(size(B)); 

for idim=1:ndim 
    rowB = idim:ndim:ndim*nelem ; 
    rowJ = (idim-1)*nelem+1: idim*nelem ; 
    J(rowJ,:) = B(rowB,:) ; 
end

j = J' ; 
j = j(:) ;