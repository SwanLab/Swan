function   XeALL= COORallELEM(ndim,nelem,nnodeE,CN,COOR)
% COORallELEM returns a matrix XeALL such that 
%  J = [J_elem1 ; J_elem2 ... ] = XeALL*BlocXiM' , that is,
% it allows one to compute the Jacobian matrix at any Gauss point by simply
% multiplying by the matrix of gradients of shape functions
XeALL = zeros(ndim*nelem,nnodeE);
for inode = 1:nnodeE
    setnodes = CN(:,inode) ; % inode-th node of all elements
    for idim = 1:ndim
        setDOFS = idim:ndim:nelem*ndim ;   % Global DOFs   
        XeALL(setDOFS,inode) = COOR(setnodes,idim) ;
    end
end