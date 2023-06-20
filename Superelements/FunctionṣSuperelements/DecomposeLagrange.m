function [lambda,initialDof] = DecomposeLagrange(GlobalU,boundaries,subC,initialDof,type)

    if type == "Two"
        param = 1;
        lambda = cell(1,boundaries); 
    elseif type == "Three"
        param = 2;
        lambda = cell(1,2*boundaries - 1); 
    end

    for i=1:boundaries
        for j=1:param
            if (i==boundaries && j==2)
                % This boundary is not connected to anything.
            else
            DofsBoundLagrange = (1:1:size(subC{i,j},2)) + initialDof;
            subLambda = GlobalU(DofsBoundLagrange);
            subLambda = [subLambda(1:2:end-1), subLambda(2:2:end)];
            lambda(param*i+j-param) = {subLambda};
            initialDof = max(DofsBoundLagrange);
            end
        end
    end
end