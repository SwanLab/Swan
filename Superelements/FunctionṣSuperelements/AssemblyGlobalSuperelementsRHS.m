function [globalRHS] = AssemblyGlobalSuperelementsRHS(RHSinputs)
    subK = RHSinputs.subK;
    force = RHSinputs.force;
    dofsBoundary = RHSinputs.dofsBoundary;
    globalLHS = RHSinputs.globalLHS;
    
    sumDofsDomain = 0;
    for i=1:(length(subK)-1)
        sumDofsDomain = sumDofsDomain + size(subK{i},1);
    end
    totalDofsDomain = sumDofsDomain + size(subK{end},1);
    nDofsLeft = size(dofsBoundary(:,1),1);
    
    dofsNeumann = (sumDofsDomain+nDofsLeft+force.dim):2:totalDofsDomain;
    nDofsNeumann = length(dofsNeumann);

    globalRHS = sparse(size(globalLHS,1),1);
    globalRHS(dofsNeumann) = force.magnitude/nDofsNeumann;
end



