function [globalRHS] = AssemblyGlobalRHS(RHSinputs) 
        subBoundMesh = RHSinputs.subBoundMesh;
        force = RHSinputs.force;
        subK = RHSinputs.subK;
        globalLHS = RHSinputs.globalLHS;

        sumDofsSubdomain = 0;
        for i=1:(length(subK)-1)
            sumDofsSubdomain = sumDofsSubdomain + size(subK{i},1);
        end
        ndim = subBoundMesh(end).mesh.ndim;
        force.dim = 2; 
    
        idxNeumann = [subBoundMesh(end,2).globalConnec(:,1); subBoundMesh(end,2).globalConnec(end,2)];
        dofsNeumann = ndim*idxNeumann+(-ndim+force.dim) + sumDofsSubdomain;
        nDofsNeumann = length(dofsNeumann);

        globalRHS = sparse(size(globalLHS,1),1);
        globalRHS(dofsNeumann) = force.magnitude/nDofsNeumann;
end