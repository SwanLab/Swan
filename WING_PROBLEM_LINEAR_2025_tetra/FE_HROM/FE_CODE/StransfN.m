function Ne = StransfN(NeSCL,ndim)
nnodeE = size(NeSCL,2); 
Ne = zeros(ndim, ndim*nnodeE); 
for inode = 1:nnodeE
    ini = (inode-1)*ndim + 1; fin = inode*ndim ; 
    Ne(:,ini:fin) = NeSCL(inode)*eye(ndim) ; 
end