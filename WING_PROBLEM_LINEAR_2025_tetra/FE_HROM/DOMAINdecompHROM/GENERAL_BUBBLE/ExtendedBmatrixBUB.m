function Bmat = ExtendedBmatrixBUB(BmatB,BmatBUB,ndim,nBUB)


if ndim ==nBUB
    % Nothing is to be done
    Bmat = [BmatB,BmatBUB] ;
elseif nBUB >  ndim
    % Number of bubble modes greater than the number spatial dimensions
    % Boundary interscale B-matrix is to be expanded with zeros
    nnodeE_b = size(BmatB,2)/ndim ;  % Number of boundary nodes
    nDOFStotB = nnodeE_b*nBUB;   % Total number of DOFs, including the "ghost" dofs
    nrows = size(BmatB,1) ; % Number of ECM points times the number of strain components
    BmatB_exp = zeros(nrows,nDOFStotB) ;
    for idim = 1:nnodeE_b
        iiniNEW = (idim-1)*nBUB + 1;
        ifinNEW = (idim-1)*nBUB + ndim ;
        iiniOLD = (idim-1)*ndim + 1;
        ifinOLD = (idim-1)*ndim + ndim ;
        BmatB_exp(:,iiniNEW:ifinNEW) = BmatB(:,iiniOLD:ifinOLD) ;
    end
    Bmat = [BmatB_exp,BmatBUB] ;
else
    % Number of bubble modes less than number of spatial dimensions
    %  nnodeE_b = size(BmatB,2)/ndim ;  % Number of boundary nodes
    diffZEROS = ndim-nBUB ;
    nrows = size(BmatB,1) ;
    cZ = zeros(nrows,diffZEROS) ;
    Bmat = [BmatB,BmatBUB,cZ] ;
end