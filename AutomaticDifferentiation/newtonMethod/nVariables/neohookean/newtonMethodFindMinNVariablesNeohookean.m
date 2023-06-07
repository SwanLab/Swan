function [plotU,plotVal,plotGrad,plotGrad2,iterations] = newtonMethodFindMinNVariablesNeohookean(u)

iterations = 10;

[~,numElem] = size(u);

plotU = zeros(iterations,numElem);
plotVal = zeros(iterations,1);
plotGrad = zeros(iterations,numElem);
plotGrad2 = zeros(iterations,numElem);

iterations = 0;

grad = 1 + zeros(1,numElem);

while abs(sum(grad(1,2:end))) > 10^(-6) && iterations < 10^4
    iterations = iterations + 1;
    [val, grad, grad2] = iterativeADNewtonNVariablesNeohookean(u);
    h = grad2.^(-1);
    h(isinf(h)) = 1;
    h = transpose(diag(h));
    u = u - h .* grad;
    u(1)=0;

    plotU(iterations,1:numElem) = u;
    plotVal(iterations) = val;
    plotGrad(iterations,1:numElem) = grad;
    plotGrad2(iterations,1:numElem) = transpose(diag(grad2));
end