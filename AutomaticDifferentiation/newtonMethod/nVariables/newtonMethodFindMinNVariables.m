function [plotU,plotVal,plotGrad,plotGrad2,iterations] = newtonMethodFindMinNVariables(u)

iterations = 41;

[~,numElem] = size(u);

plotU = zeros(iterations,numElem);
plotVal = zeros(iterations,1);
plotGrad = zeros(iterations,numElem);
plotGrad2 = zeros(iterations,numElem);

iterations = 0;

grad = 1 + zeros(1,numElem);

while abs(sum(abs(grad(1,3:end)))) > 10^(-6) && iterations < 10^5
    iterations = iterations + 1;
    [val, grad, grad2] = iterativeADNewtonNVariables(u);
    h = grad2.^(-1);
    %h(isinf(h)) = 0;
    h = transpose(diag(h));
    u = u - h .* grad;
    u(1) =  0;

    plotU(iterations,1:numElem) = u;
    plotVal(iterations) = val;
    plotGrad(iterations,1:numElem) = grad;
    plotGrad2(iterations,1:numElem) = transpose(diag(grad2));
end


%%PLOTS

% %figure(1); plot(plotU(:,1)); xlabel("Num. of iterations"); ylabel("Displacement"); grid; axis([0 iterations min(plotU(:,1)) max(plotU(:,1))]);
% figure(2); plot(plotU(:,2)); xlabel("Num. of iterations"); ylabel("Displacement"); grid; axis([0 iterations min(plotU(:,2)) max(plotU(:,2))]);
% figure(3); plot(plotU(:,3)); xlabel("Num. of iterations"); ylabel("Displacement"); grid; axis([0 iterations min(plotU(:,3)) max(plotU(:,3))]);
% %figure(4); plot(plotGrad(:,1)); xlabel("Num. of iterations"); ylabel("First Gradient"); grid; axis([0 iterations min(plotGrad(:,1)) max(plotGrad(:,1))]);
% figure(5); plot(plotGrad(:,2)); xlabel("Num. of iterations"); ylabel("First Gradient"); grid; axis([0 iterations min(plotGrad(:,2)) max(plotGrad(:,2))]);
% figure(6); plot(plotGrad(:,3)); xlabel("Num. of iterations"); ylabel("First Gradient"); grid; axis([0 iterations min(plotGrad(:,3)) max(plotGrad(:,3))]);
% %figure(7); plot(plotGrad2(:,1)); xlabel("Num. of iterations"); ylabel("Second Gradient"); grid; axis([0 iterations min(plotGrad2(:,1)) max(plotGrad2(:,1))]);
% %figure(8); plot(plotGrad2(:,2)); xlabel("Num. of iterations"); ylabel("Second Gradient"); grid; axis([0 iterations min(plotGrad2(:,2)) max(plotGrad2(:,2))]);
% %figure(9); plot(plotGrad2(:,3)); xlabel("Num. of iterations"); ylabel("Second Gradient"); grid; axis([0 iterations min(plotGrad2(:,3)) max(plotGrad2(:,3))]);
% figure(10); plot(plotVal); xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations min(plotVal) max(plotVal)]);