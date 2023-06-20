function [U,initialDof] = DecomposeSuperelementsDisplacement(GlobalU,subdomains,dofsBoundary)
    U = cell(subdomains,2);
    initialDof = 0;
    for i=1:subdomains
        for j=1:2
            if j==1
                nDofsBoundary = length(dofsBoundary(:,1));
            elseif j==2
                nDofsBoundary = length(dofsBoundary(:,2));
            end
            DofsBoundary = (1:1:nDofsBoundary) + initialDof;
            subU = GlobalU(DofsBoundary);
            subU = [subU(1:2:end-1), subU(2:2:end)];
            U(i,j) = {subU};
            initialDof = max(DofsBoundary);
        end
    end
end