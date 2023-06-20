
function [globalC] = AssemblyGlobalLagrange(ACinputs)
    type = ACinputs.type;
    subBoundMesh = ACinputs.subBoundMesh;
    subC = ACinputs.subC;
    subK = ACinputs.subK;
    dofU = size(ACinputs.globalK,1);

    boundaries = size(subC,1);
    ndim = subBoundMesh(1).mesh.ndim;

    % Start of the assembly
    sumDofSubdomains = 0;
    dofSubLagrange =  size(subC{1,1},2);
    columnC = sparse(dofU,dofSubLagrange);

    idx = [subBoundMesh(1,1).globalConnec(:,1); subBoundMesh(1,1).globalConnec(end,2)];
    dofsBoundary = computeDofsBoundary(idx,ndim,sumDofSubdomains);

    columnC(dofsBoundary,:) = -subC{1,1};
    globalC = columnC;

    for i=2:boundaries
        for j=1:2
            idxSubdomain = i+j-2;
            idxBoundary  = 3-j;

            dofSubdomain = size(subK{i},1);
            dofSubLagrange = size(subC{i,1},2);
            if j==1 || (j==2 && type == "Three")
            columnC = sparse(dofU,dofSubLagrange);
            end
            idx = [subBoundMesh(idxSubdomain,idxBoundary).globalConnec(:,1); subBoundMesh(idxSubdomain,idxBoundary).globalConnec(end,2)];
            dofsBoundary = computeDofsBoundary(idx,ndim,sumDofSubdomains);

            if j==1
                columnC(dofsBoundary,:) = subC{idxSubdomain,idxBoundary};
            elseif j==2
                columnC(dofsBoundary,:) = -subC{idxSubdomain,idxBoundary};
            end

            if (j==1 && type == "Three") || j==2
                globalC = [globalC columnC];
            end
            sumDofSubdomains = sumDofSubdomains + dofSubdomain;
        end
        sumDofSubdomains = sumDofSubdomains - dofSubdomain;
    end
end




