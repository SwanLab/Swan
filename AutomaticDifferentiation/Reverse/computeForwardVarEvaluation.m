function [time,grad] = computeForwardVarEvaluation(values)
    tic
    n = length(values);
    x = cell(1,n);
    for i = 1:n
        grad = zeros(1,n);
        grad(i) = 1;
        x{i} = ValGradForward(values(i),grad);
    end
    f = 0;
    for i = 2:n-1
        f = f + x{i}^x{i+1}*cos(i-1)*log(i);
    end
    AD = f.double;
    val = AD(:,1);
    grad = AD(:,2:end);
    time = toc;
end