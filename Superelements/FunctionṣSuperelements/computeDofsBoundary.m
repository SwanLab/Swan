function [dofsBoundary] = computeDofsBoundary(idx,ndim,sumDofsSubdomain)
    dofsBoundary = zeros(length(idx)*ndim,1);
    for i=1:ndim
        dofsBoundary(i:2:(end-ndim+i)) = ndim*idx+(-ndim+i) + sumDofsSubdomain;
    end
end