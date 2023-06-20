function [uFull] = FullMeshReconstruction(Reconstruct)
    globalU = Reconstruct.globalU;
    superK = Reconstruct.superElemInfo.superK;
    dofsCondensed = Reconstruct.superElemInfo.dofs.condensed;
    dofsBoundary = Reconstruct.superElemInfo.dofs.boundary;
    dofsBoundary = [dofsBoundary(:,1); dofsBoundary(:,2)];
    K = Reconstruct.superElemInfo.K;

    subdomains = length(superK);
    uFull = cell(1,subdomains);
    for i=1:subdomains
        uB = globalU((i-1)*length(dofsBoundary)+1:i*length(dofsBoundary));
        uI = K.ii\(-K.ib*uB);
    
        u = zeros(length(uB)+length(uI),1);
        u(dofsCondensed) = uI;
        u(dofsBoundary) = uB;

        u = [u(1:2:end-1), u(2:2:end)];
        uFull(i) = {u};
    end
end