function [U,initialDof] = DecomposeDisplacement(GlobalU,subdomains,subK)
    U = cell(1,subdomains);
    initialDof = 0;
    for i=1:subdomains
        DofsSubdomain = (1:1:size(subK{i},1)) + initialDof;
        subU = GlobalU(DofsSubdomain);
        subU = [subU(1:2:end-1), subU(2:2:end)];
        U(i) = {subU};
        initialDof = max(DofsSubdomain);
    end
end