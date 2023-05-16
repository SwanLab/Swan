function [plotU,plotVal,plotGrad,iterations] = newtonMethodFindMin(u)
iterations = 36;
plotU = zeros(iterations,1);
plotVal = zeros(iterations,1);
plotGrad = zeros(iterations,1);
alpha = 0.1;
iterations = 0;
grad = 1;
while abs(grad) > 10^(-6)
    iterations = iterations + 1;
    [val, grad] = iterativeADNewton(u);
    u = u - alpha * grad;
    plotU(iterations) = u;
    plotVal(iterations) = val;
    plotGrad(iterations) = grad;
end

figure(1);
plot(plotU); xlabel("Num. of iterations"); ylabel("Displacement"); grid; axis([0 iterations min(plotU) max(plotU)]);
figure(2);
plot(plotVal); xlabel("Num. of iterations"); ylabel("Value"); grid; axis([0 iterations min(plotVal) max(plotVal)]);
figure(3);
plot(plotGrad); xlabel("Num. of iterations"); ylabel("Gradient"); grid; axis([0 iterations min(plotGrad) max(plotGrad)]);
figure(4);
plot(plotU,plotVal); grid; axis([min(plotU) max(plotU) min(plotVal) max(plotVal)]);