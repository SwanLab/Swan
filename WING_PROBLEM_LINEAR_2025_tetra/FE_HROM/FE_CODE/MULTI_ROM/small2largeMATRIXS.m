function M = small2largeMATRIXS(M1d,ndim)

 M = sparse(size(M1d,1)*ndim,size(M1d,2)*ndim) ;
    for idim = 1:ndim
        M(idim:ndim:end,idim:ndim:end) = M1d ;
    end
