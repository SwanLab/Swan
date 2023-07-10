function [time,grad] = computeReverseVarEvaluation(values)
    tic
    n = length(values);
    x = cell(1,n);
    for i = 1:n
        x{i} = ReverseVar(values(i));
    end
    f = 0;
    for i = 2:n-1
        f = f + x{i}^x{i+1}*cos(i-1)*log(i);
    end
    f.grad_value = 1;
    grad = cell(1,n);
    for i = 1:n
        x{i}.computeGradient;
        grad{i} = x{i}.grad_value;
    end
    grad = cell2mat(grad);
    time = toc;
end