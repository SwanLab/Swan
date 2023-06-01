function [plotU,plotVal,plotGrad,plotGrad2,iterations] = newtonMethodFindMinNVariablesPlotDisp(u)

nodePlace = zeros(1,length(u));

for i = 1:length(u)
    nodePlace(i) = i-1;
end

iterations = 917;

[~,numElem] = size(u);

plotU = zeros(iterations,numElem);
plotVal = zeros(iterations,1);
plotGrad = zeros(iterations,numElem);
plotGrad2 = zeros(iterations,numElem);

iterations = 0;

grad = 1 + zeros(1,numElem);

while abs(sum(grad(1,3:end))) > 10^(-6) && iterations < 10^5
    iterations = iterations + 1;
    [val, grad, grad2] = iterativeADNewtonNVariablesPlotDisp(u);
    h = grad2.^(-1);
    h(isinf(h)) = 0;
    h = transpose(diag(h));
    u = u - h .* grad;

    plotU(iterations,1:numElem) = u;
    plotU(iterations,1) = 0;
    plotVal(iterations) = val;
    plotGrad(iterations,1:numElem) = grad;
    plotGrad2(iterations,1:numElem) = transpose(diag(grad2));

    figure(1);
    plot(nodePlace,-0.5+zeros(1,length(u)),'o',nodePlace(1),0,'o'); grid; axis([0 nodePlace(length(u))+u(length(u)) -1 1]);
    hold on
    for i = 2:length(u)
        plot(nodePlace(i)+u(i),0,'o');
    end
    hold off


end