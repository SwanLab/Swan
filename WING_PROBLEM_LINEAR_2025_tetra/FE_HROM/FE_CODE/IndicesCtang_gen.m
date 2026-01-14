function   [i,j] = IndicesCtang_gen(m,p,nrows)

% See ConvertBlockDiag_general.m

if nargin == 0
    p = 2;
    m = 6 ;
end


ngaus = m/nrows ;
n = ngaus*p ; 
nzmax = m*p ;
ij = zeros(p*nrows,2) ;
for ielem = 1:nrows
    for jelem = 1:p
        iniELEM = (ielem-1)*p +jelem;
        ij(iniELEM,1) = ielem ;
        ij(iniELEM,2) = jelem ;
    end
end

ij = repmat(ij,ngaus,1);
indSUM = 0:p:(ngaus-1)*p ;
indSUM = repmat(indSUM,p^2,1);
indSUM = reshape(indSUM,size(indSUM,1)*size(indSUM,2),1) ;
ij = bsxfun(@plus,ij,indSUM) ;
i = ij(:,1);
j = ij(:,2);