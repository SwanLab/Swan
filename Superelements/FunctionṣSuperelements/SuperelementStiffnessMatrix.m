function [SuperElemInfo] = SuperelementStiffnessMatrix(CondensInputs)
    subK = CondensInputs.subK;
    subBoundMesh = CondensInputs.subBoundMesh;

    subdomains = length(subK);
    ndim = subBoundMesh(1,1).mesh.ndim;
    dofsBoundary = [];
    for j=1:2
        idx = [subBoundMesh(1,j).globalConnec(:,1); subBoundMesh(1,j).globalConnec(end,2)];
        dofsBoundary = [dofsBoundary, computeDofsBoundary(idx,ndim,0)];
    end
    dofsBoundSuper = [dofsBoundary(:,1); dofsBoundary(:,2)];
    totalDofs = 1:1:size(subK{1},1);
    dofsCondensed = setdiff(totalDofs,dofsBoundSuper)';

    dofs.condensed = dofsCondensed;
    dofs.boundary = dofsBoundary;

    K.ii = subK{1}(dofsCondensed,dofsCondensed);
    K.ib = subK{1}(dofsCondensed,dofsBoundSuper);
    K.bi = subK{1}(dofsBoundSuper,dofsCondensed);
    K.bb = subK{1}(dofsBoundSuper,dofsBoundSuper);

    superK = K.bb - K.bi*(K.ii\K.ib);
    for k=1:subdomains 
        subK{k} = superK;
    end

    SuperElemInfo.K = K;
    SuperElemInfo.superK = subK;
    SuperElemInfo.dofs = dofs;
end
